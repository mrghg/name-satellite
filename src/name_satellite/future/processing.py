"""Future-method processing: NAME output -> column footprints on the nested grids, one file a month.

Where this sits (docs/future/scoping.md; the file layout: docs/future/file-format.md)
------------------------------------------------------------------------------------
    ``future.setup`` -> NAME on LOTUS2 -> **future.processing** -> the inversion

Per day (``process_day``), from the run folder written by ``future.setup``:

1. Inputs, matched by ID: the observation file (sounding IDs, weights), ``release_groups.csv``
   (sounding -> release group), ``release_plan.csv`` (releases, levels), ``grids.csv`` (each
   group's boxes and NAME grid names), the fields files (one per output grid) and the patched
   NAME's fate records (``Particles_C1.txt``; repeats removed, X20).
2. Per release (level j of group g): its footprint on the group's 0.1° box, 0.5° box and the
   global 2.5° grid, (mol/mol)/(mol m-2 s-1): field (ppm s) x 1e-6 x cell area / release
   duration (``footprints``); and where its particles ended: a histogram on 2.5° x 1 km bins
   (0-20 km agl) plus one bin for particles that left through the 20 km top, as **fractions of
   the particles released** (counted from the fate records), so each release's histogram sums
   to 1 (checked).
3. Per sounding s: every array = Σ_j w_sj × (its group's level-j array), w = pressure weight ×
   averaging kernel over the modelled levels, release level j = retrieval level j (as the main
   processing, D23/D26: one column footprint per sounding). The end-position histogram then sums
   to Σ_j w_sj (checked).
4. The 0.1° and 0.5° boxes go into fixed-size blocks (``future.grids``: 125 x 125, 105 x 105),
   zero-padded, with their offset in the global grid of that level and their true shape.

``process_month`` assembles the days, sorted by time and sounding ID, and ``write_file`` writes
float32 arrays rounded to 10 significant bits (BitRound) and compressed (zlib 4, shuffle), one
sounding per chunk, as the main processing does.

Command line::

    uv run python -m name_satellite.future.processing MONTH_FOLDER OUT_FOLDER [--days=YYYYMMDD,...]
        [--allow-test]
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from name_satellite.footprints import footprint_from_mixing_ratio, grid_cell_area_m2
from name_satellite.future.grids import (
    COARSE,
    COARSE_NX,
    COARSE_NY,
    FINE,
    MEDIUM,
    Box,
    box_offset,
    box_shape,
    global_centres,
    offset_from_centres,
)
from name_satellite.future.setup import MAX_AGE_HOURS, SYNC_MINUTES
from name_satellite.name_output import (
    drop_repeated_fate_records,
    read_fields,
    read_name_fate_file,
)
from name_satellite.processing import software_commit
from name_satellite.releases import MODELLED_LEVELS

END_DZ_M = 1000.0
END_TOP_M = 20000.0
N_END_HEIGHTS = int(END_TOP_M / END_DZ_M) + 1  # 20 bins of 1 km, then "left through the top"
LEFT_DOMAIN = 8  # fate-code bit (patched NAME, X20): the particle left the domain (here: the top)
BITROUND_BITS = 10
ROUNDED = ("srr_fine", "srr_medium", "srr_coarse", "end_positions")


# ---------------------------------------------------------------------------------------------
# Reading a run
# ---------------------------------------------------------------------------------------------


def read_run(day_folder, allow_test=False):
    """Everything ``process_day`` needs from one day's run folder (see the module docstring).

    Raises:
        ValueError for a test run (unless ``allow_test``) or a missing fields file.
    """
    day_folder = Path(day_folder)
    info = json.loads((day_folder / "day.json").read_text())
    if info["options"]["test"] and not allow_test:
        raise ValueError(f"{day_folder} is a test run (1 h of travel); use --allow-test")
    day = pd.Timestamp(info["day"])
    with xr.open_dataset(day_folder / f"obs_{info['site']}_{day:%Y%m%d}.nc") as f:
        obs = f.load()
    plan = pd.read_csv(day_folder / "release_plan.csv", dtype={"release_group": str},
                       parse_dates=["release_start", "release_end"])  # fmt: skip
    members = pd.read_csv(day_folder / "release_groups.csv", dtype=str)
    boxes = pd.read_csv(day_folder / "grids.csv", dtype={"release_group": str})
    fields = {}
    for grid in [*boxes["fine_grid"], *boxes["medium_grid"], "HGridC"]:
        files = list((day_folder / "output").glob(f"Fields_{grid}_C1_T1_*.txt"))
        if len(files) != 1:
            raise ValueError(f"{day_folder}: expected one fields file for {grid}, found {files}")
        fields[grid] = read_fields(files[0])
    particles, z_units = read_name_fate_file(day_folder / "output" / "Particles_C1.txt")
    if z_units != "m agl":
        raise ValueError(f"Fate heights in {z_units}, expected m agl")
    # Repeated exit records lie at most one sync step apart (X20).
    particles = drop_repeated_fate_records(particles, max_gap_seconds=60.0 * SYNC_MINUTES)
    return {"info": info, "obs": obs, "plan": plan, "members": members, "boxes": boxes,
            "fields": fields, "particles": particles}  # fmt: skip


# ---------------------------------------------------------------------------------------------
# Per release
# ---------------------------------------------------------------------------------------------


def level_footprints(field, sources, release_seconds):
    """Footprints (level, lat, lon) of the given sources on one grid, (mol/mol)/(mol m-2 s-1)."""
    area = grid_cell_area_m2(field["lat"].values, field["lon"].values)
    return footprint_from_mixing_ratio(field.sel(source=sources).values, area, release_seconds)


def end_position_fractions(particles):
    """Histogram of one release's end positions as fractions of its particles.

    Args:
        particles: fate records of one source (one per particle), with lon, lat, z (m agl),
            fate_code.

    Returns:
        Array (N_END_HEIGHTS, COARSE_NY, COARSE_NX) summing to 1: heights 0-20 km in 1 km bins,
        the last level = left through the top (fate code with bit 8); cells of the global 2.5°
        grid.
    """
    n = len(particles)
    if n == 0:
        raise ValueError("A release has no fate records")
    i = np.floor((particles["lon"].to_numpy() + 180.0) / COARSE.cell_deg).astype(int) % COARSE_NX
    j = np.clip(np.floor((particles["lat"].to_numpy() + 90.0) / COARSE.cell_deg).astype(int), 0,
                COARSE_NY - 1)  # fmt: skip
    k = np.clip(np.floor(particles["z"].to_numpy() / END_DZ_M).astype(int), 0, N_END_HEIGHTS - 2)
    left_top = (particles["fate_code"].to_numpy() & LEFT_DOMAIN) != 0
    k[left_top] = N_END_HEIGHTS - 1
    hist = np.zeros((N_END_HEIGHTS, COARSE_NY, COARSE_NX))
    np.add.at(hist, (k, j, i), 1.0 / n)
    return hist


def check_box(field, box, level):
    """The NAME grid of a box is the snapped box (shape and global offset)."""
    lat, lon = field["lat"].values, field["lon"].values
    if (len(lat), len(lon)) != box_shape(box, level.cell_deg):
        raise ValueError(f"{level.name} grid shape {(len(lat), len(lon))} differs from its box")
    if offset_from_centres(lat, lon, level.cell_deg) != box_offset(box, level.cell_deg):
        raise ValueError(f"{level.name} grid is not at its box's place in the global grid")


def in_block(values, level):
    """Values (..., ny, nx) zero-padded into the level's fixed block (..., n, n)."""
    n = level.block_cells
    out = np.zeros((*values.shape[:-2], n, n))
    out[..., : values.shape[-2], : values.shape[-1]] = values
    return out


# ---------------------------------------------------------------------------------------------
# One day
# ---------------------------------------------------------------------------------------------


def process_day(day_folder, allow_test=False):
    """Column footprints and end positions of every sounding of one day (xarray.Dataset).

    Raises:
        ValueError if IDs disagree, a group lacks levels 1..n, a grid is not at its box, or the
        particle accounting does not add up.
    """
    run = read_run(day_folder, allow_test)
    obs, plan, members, boxes = run["obs"], run["plan"], run["members"], run["boxes"]
    fields, particles = run["fields"], run["particles"]
    product = run["info"]["product"]
    n_levels = MODELLED_LEVELS[product]
    durations = (plan["release_start"] - plan["release_end"]).dt.total_seconds().unique()
    if len(durations) != 1:
        raise ValueError(f"Releases of different durations: {durations}")
    release_seconds = float(durations[0])
    max_age_hours = 1 if run["info"]["options"]["test"] else MAX_AGE_HOURS

    sounding_ids = [str(s) for s in obs["sounding_id"].values]
    group_of = members.set_index("sounding_id")["release_group"]
    if sorted(group_of.index) != sorted(sounding_ids):
        raise ValueError("release_groups.csv and the observation file have different soundings")
    if set(group_of) != set(boxes["release_group"]) or set(group_of) != set(plan["release_group"]):
        raise ValueError("Release groups differ between release_groups.csv, plan and grids.csv")
    weights = (obs["pressure_weights"].values[:, :n_levels]
               * obs["xch4_averaging_kernel"].values[:, :n_levels]).astype(float)  # fmt: skip
    if not np.isfinite(weights).all():
        raise ValueError("Missing pressure weights or averaging kernels in the modelled levels")
    particles["source_number"] = particles["source_number"].astype(int)
    by_source = dict(tuple(particles.groupby("source_number")))

    per_group = {}
    for row in boxes.to_dict("records"):
        releases = plan[plan["release_group"] == row["release_group"]].sort_values("level_number")
        if list(releases["level_number"]) != list(range(1, n_levels + 1)):
            raise ValueError(f"Group {row['release_group']} lacks levels 1..{n_levels}")
        sources = list(releases["source_name"])
        fine_box = Box(row["fine_west"], row["fine_east"], row["fine_south"], row["fine_north"])
        medium_box = Box(row["medium_west"], row["medium_east"], row["medium_south"],
                         row["medium_north"])  # fmt: skip
        check_box(fields[row["fine_grid"]], fine_box, FINE)
        check_box(fields[row["medium_grid"]], medium_box, MEDIUM)
        ends = []
        for source in sources:
            number = int(source.removeprefix("SourceID_"))
            if number not in by_source:
                raise ValueError(f"No fate records for {source}")
            ends.append(end_position_fractions(by_source[number]))
        ends = np.array(ends)
        if not np.allclose(ends.sum(axis=(1, 2, 3)), 1.0, rtol=0, atol=1e-9):
            raise ValueError(f"End positions of group {row['release_group']} do not sum to 1")
        first = releases.iloc[0]
        per_group[row["release_group"]] = {
            "fine": level_footprints(fields[row["fine_grid"]], sources, release_seconds),
            "medium": level_footprints(fields[row["medium_grid"]], sources, release_seconds),
            "coarse": level_footprints(fields["HGridC"], sources, release_seconds),
            "ends": ends,
            "fine_box": fine_box,
            "medium_box": medium_box,
            "release": (
                float(first["lon"]),
                float(first["lat"]),
                first["release_start"],
                first["release_end"],
            ),
            "top_age": [by_source[int(s.removeprefix("SourceID_"))] for s in sources],
        }

    data = {name: [] for name in ("srr_fine", "srr_medium", "srr_coarse", "end_positions")}
    extra = {name: [] for name in ("fine_offset", "fine_shape", "medium_offset", "medium_shape",
                                   "release_lon", "release_lat", "release_time", "end_time",
                                   "top_exit_mean_age")}  # fmt: skip
    for s, sounding in enumerate(sounding_ids):
        g = per_group[group_of[sounding]]
        w = weights[s]
        data["srr_fine"].append(in_block(np.tensordot(w, g["fine"], axes=1), FINE))
        data["srr_medium"].append(in_block(np.tensordot(w, g["medium"], axes=1), MEDIUM))
        data["srr_coarse"].append(np.tensordot(w, g["coarse"], axes=1))
        ends = np.tensordot(w, g["ends"], axes=1)
        if not np.isclose(ends.sum(), w.sum(), rtol=1e-9, atol=0):
            raise ValueError(f"End positions of sounding {sounding} do not add up to its weight")
        data["end_positions"].append(ends)
        extra["fine_offset"].append(box_offset(g["fine_box"], FINE.cell_deg))
        extra["fine_shape"].append(box_shape(g["fine_box"], FINE.cell_deg))
        extra["medium_offset"].append(box_offset(g["medium_box"], MEDIUM.cell_deg))
        extra["medium_shape"].append(box_shape(g["medium_box"], MEDIUM.cell_deg))
        lon, lat, start, end = g["release"]
        extra["release_lon"].append(lon)
        extra["release_lat"].append(lat)
        extra["release_time"].append(start)
        extra["end_time"].append(end - pd.Timedelta(hours=max_age_hours))
        # Mean travel time of the particles that left through the top, weighted as the column.
        top_weight = top_age = 0.0
        for level, records in enumerate(g["top_age"]):
            top = records[(records["fate_code"].to_numpy() & LEFT_DOMAIN) != 0]
            top_weight += w[level] * len(top) / len(records)
            top_age += w[level] * top["age_hours"].sum() / len(records)
        extra["top_exit_mean_age"].append(top_age / top_weight if top_weight > 0 else np.nan)

    fine_lat, fine_lon = global_centres(FINE.cell_deg)
    medium_lat, medium_lon = global_centres(MEDIUM.cell_deg)
    coarse_lat, coarse_lon = global_centres(COARSE.cell_deg)
    heights = np.append(np.arange(END_DZ_M / 2, END_TOP_M, END_DZ_M), np.nan)
    s_dim = "sounding_id"
    return xr.Dataset(
        {
            "srr_fine": ((s_dim, "fine_y", "fine_x"), np.array(data["srr_fine"])),
            "fine_offset": ((s_dim, "yx"), np.array(extra["fine_offset"], dtype="int32")),
            "fine_shape": ((s_dim, "yx"), np.array(extra["fine_shape"], dtype="int32")),
            "srr_medium": ((s_dim, "medium_y", "medium_x"), np.array(data["srr_medium"])),
            "medium_offset": ((s_dim, "yx"), np.array(extra["medium_offset"], dtype="int32")),
            "medium_shape": ((s_dim, "yx"), np.array(extra["medium_shape"], dtype="int32")),
            "srr_coarse": ((s_dim, "coarse_lat", "coarse_lon"), np.array(data["srr_coarse"])),
            "end_positions": (
                (s_dim, "end_height", "coarse_lat", "coarse_lon"),
                np.array(data["end_positions"]),
            ),
            "column_weight": ((s_dim, "release_level"), weights),
            "release_group": (s_dim, np.array([group_of[s] for s in sounding_ids])),
            "time": (s_dim, obs["time"].values),
            "lat": (s_dim, obs["lat"].values.astype(float)),
            "lon": (s_dim, obs["lon"].values.astype(float)),
            "release_lon": (s_dim, np.array(extra["release_lon"])),
            "release_lat": (s_dim, np.array(extra["release_lat"])),
            "release_time": (s_dim, pd.to_datetime(extra["release_time"]).values),
            "end_time": (s_dim, pd.to_datetime(extra["end_time"]).values),
            "top_exit_mean_age": (s_dim, np.array(extra["top_exit_mean_age"])),
        },
        coords={
            s_dim: np.array(sounding_ids),
            "release_level": np.arange(1, n_levels + 1),
            "fine_lat_global": ("fine_lat_global", fine_lat),
            "fine_lon_global": ("fine_lon_global", fine_lon),
            "medium_lat_global": ("medium_lat_global", medium_lat),
            "medium_lon_global": ("medium_lon_global", medium_lon),
            "coarse_lat": coarse_lat,
            "coarse_lon": coarse_lon,
            "end_height": heights,
        },
        attrs={"max_age_hours": max_age_hours, "options": json.dumps(run["info"]["options"])},
    )


# ---------------------------------------------------------------------------------------------
# A month, the file
# ---------------------------------------------------------------------------------------------


def process_month(month_folder, days=None, allow_test=False):
    """Every day folder (YYYYMMDD) of the month, or only ``days``, assembled and sorted."""
    folders = sorted(p for p in Path(month_folder).iterdir()
                     if p.is_dir() and p.name.isdigit() and len(p.name) == 8)  # fmt: skip
    if days is not None:
        folders = [p for p in folders if p.name in days]
    if not folders:
        raise ValueError(f"No day folders to process in {month_folder}")
    parts = [process_day(f, allow_test) for f in folders]
    ds = xr.concat(parts, dim="sounding_id", data_vars="minimal", coords="minimal",
                   compat="override", combine_attrs="override")  # fmt: skip
    order = np.lexsort((ds["sounding_id"].values.astype(str), ds["time"].values))
    return ds.isel(sounding_id=order)


ATTRIBUTES = {
    "srr_fine": ("Column footprint on the sounding's 0.1° box (fixed block, zero-padded; see "
                 "fine_offset, fine_shape)", "(mol/mol) / (mol m-2 s-1)"),
    "srr_medium": ("Column footprint on the sounding's 0.5° box (fixed block, zero-padded; see "
                   "medium_offset, medium_shape)", "(mol/mol) / (mol m-2 s-1)"),
    "srr_coarse": ("Column footprint on the global 2.5° grid", "(mol/mol) / (mol m-2 s-1)"),
    "fine_offset": ("Global 0.1° index (lat, lon) of the block's first cell", "1"),
    "fine_shape": ("Cells (lat, lon) of the block that belong to the box", "1"),
    "medium_offset": ("Global 0.5° index (lat, lon) of the block's first cell", "1"),
    "medium_shape": ("Cells (lat, lon) of the block that belong to the box", "1"),
    "end_positions": ("Column-weighted fraction of particles ending in each 2.5° x 1 km bin "
                      "(last height: left through the 20 km top); sums to the column weight",
                      "1"),
    "column_weight": ("Pressure weight x averaging kernel per modelled level", "1"),
    "end_time": ("Time at which particles reached the maximum travel time", ""),
    "top_exit_mean_age": ("Mean travel time of particles leaving through the top", "hours"),
}  # fmt: skip


def write_file(ds, path):
    """Write the month file (float32, BitRound 10 bits, zlib 4 + shuffle, one sounding a chunk)."""
    ds = ds.copy()
    for name, (long_name, units) in ATTRIBUTES.items():
        ds[name].attrs.update(long_name=long_name, **({"units": units} if units else {}))
    ds.attrs.update(
        title="NAME column footprints on nested grids (future method; draft)",
        layout="docs/future/file-format.md (name-satellite, branch future-method)",
        software="name-satellite",
        software_commit=software_commit(),
        created=pd.Timestamp.now(tz="UTC").isoformat(timespec="seconds"),
    )
    encoding = {}
    for name in ROUNDED:
        encoding[name] = {"dtype": "float32", "zlib": True, "complevel": 4, "shuffle": True,
                          "chunksizes": (1, *ds[name].shape[1:]), "quantize_mode": "BitRound",
                          "significant_digits": BITROUND_BITS}  # fmt: skip
    for name in ds.data_vars:
        if name not in encoding and ds[name].dtype.kind == "f":
            encoding[name] = {"dtype": "float64", "zlib": True}
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(path, engine="netcdf4", format="NETCDF4", encoding=encoding)
    return path


def sounding_footprint(ds, sounding_id, level):
    """One sounding's footprint on one level as a DataArray with lat/lon coordinates.

    ``level``: "fine", "medium" or "coarse". The fine and medium blocks are cut to their box.
    This is how an inversion places a sounding's blocks on the global grids.
    """
    one = ds.sel(sounding_id=sounding_id)
    if level == "coarse":
        return one["srr_coarse"].rename(coarse_lat="lat", coarse_lon="lon")
    (j0, i0), (ny, nx) = one[f"{level}_offset"].values, one[f"{level}_shape"].values
    block = one[f"srr_{level}"].values[:ny, :nx]
    lat = ds[f"{level}_lat_global"].values[j0 : j0 + ny]
    lon = ds[f"{level}_lon_global"].values[i0 : i0 + nx]
    return xr.DataArray(block, dims=("lat", "lon"), coords={"lat": lat, "lon": lon})


if __name__ == "__main__":
    arguments = [a for a in sys.argv[1:] if not a.startswith("--")]
    chosen = [a.removeprefix("--days=") for a in sys.argv[1:] if a.startswith("--days=")]
    month_folder, out_folder = Path(arguments[0]), Path(arguments[1])
    dataset = process_month(month_folder, chosen[0].split(",") if chosen else None,
                            allow_test="--allow-test" in sys.argv[1:])  # fmt: skip
    path = write_file(dataset, out_folder / f"{month_folder.name}.nc")
    print(f"Wrote {path}: {dataset.sizes['sounding_id']} soundings")
