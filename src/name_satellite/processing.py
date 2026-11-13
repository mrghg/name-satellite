"""Step 6 for new runs: one day's NAME output -> column footprints in the new file format.

Where this sits (see docs/workflow.md)
--------------------------------------
    4. NAME runs (``prepare.py``)  ->  5. read NAME output  ->  **6. processing (this module)**
    ->  footprint file read by the inversion (``docs/footprint-file-spec.md``, D23)

In (one day's run folder, written by ``prepare.py`` and NAME):
    obs_<site>_<YYYYMMDD>.nc   the soundings (keyed by sounding ID, D18)
    release_groups.csv         each sounding's release group (D20)
    release_plan.csv           one set of releases (NAME sources) per group, levels 1..16
    output/Fields_grid1*.txt   time-integrated mixing ratio of every source (ppm s)
    output/Particles_C1.txt    one fate record per particle (ACRG-patched NAME, X20)
    output/ReleaseMet_*.txt    met at each group's release point (X21)

Out: an xarray Dataset with one column footprint per sounding (``process_day``), or a month of them
(``process_month``), written by ``write_footprint_file`` as one netCDF file per month, named by
``footprint_file_name``, e.g. ``GOSAT1v10-BRAZIL-column_NAME_UMG_SOUTHAMERICA_inert_201301.nc``.
Layout (``docs/footprint-file-spec.md``): dimensions ``sounding_id`` (one entry per sounding,
sorted by time then ID), ``footprint`` (reached through ``footprint_index``), ``latitude``,
``longitude``, ``height`` (20 bins of 1 km, m above the model ground), ``release_level`` (1-16);
main variables ``srr`` ((mol/mol)/(mol/m2/s)), ``particle_locations_n/e/s/w/top``,
``mean_age_particles_*`` (hours), ``particle_fraction_inside``, the release-point met, release
layers (hPa) and ``column_weight``; provenance in the global attributes. Look at it with
``ncdump -h`` or ``xarray.open_dataset``; per footprint, the particle fractions sum to
Σ ``column_weight`` (checked here before writing).

Command line::

    uv run python -m name_satellite.processing MONTH_FOLDER OUT_FOLDER [--days=YYYYMMDD,...]

MONTH_FOLDER is ``prepare``'s month folder (every ``YYYYMMDD`` day folder in it is processed,
unless ``--days`` names some); the file is written into OUT_FOLDER and must not exist yet.
``--test`` runs are refused. A single day is fine on the sci node; a full month is heavy
processing, for LOTUS2 (docs/jasmin.md).

The steps, for each day:

1. **Per-release footprints** (6a): each source's field (ppm s) -> footprint
   ((mol/mol)/(mol/m2/s)) with ``footprints.footprint_from_mixing_ratio``.
2. **Particle fates** (6b, D13): every particle released is counted exactly once, in one fate:
   it left through the north, east, south or west face, left through the top, or was still
   inside at the maximum travel time. ``classify_fates`` says which, from its last position.
   Each face gives a histogram over (height bin, cell along the edge); the top a histogram over
   the output grid; "inside" only a count. All are divided by the particles released (counted
   from the fate records: NAME does not always release rate × duration, name-behaviour.md), so
   each release's fates sum to exactly 1. Heights are m above ground in the site files' 20 bins
   of 1 km (D24); travel times are kept as the mean age per bin.
3. **Columns** (6c): sounding s in release group g gets, for each quantity,
   Σ_j w_sj × (value of group g's release at level j),  w_sj = pressure weight × averaging
   kernel of sounding s at retrieval level j, over the levels NAME models (D21). Release level j
   represents retrieval level j by construction (the release layers are built from the
   retrieval levels, D19/D20), so they are matched by number, not by the legacy nearest-pressure
   rule (``column.nearest_release_levels``). Mean ages are combined as averages weighted by
   w_sj × (fraction of level j's particles in the bin).
4. **Met** at the group's lowest release (temperature, pressure, wind) and release point
   (boundary-layer depth, model surface pressure), from the main run (X21).

Every footprint here is exact for its sounding: ``footprint_index`` is the identity (D23).
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from name_satellite.domains import SITES
from name_satellite.footprints import footprint_from_mixing_ratio, grid_cell_area_m2
from name_satellite.name_output import (
    drop_repeated_fate_records,
    read_fields_folder,
    read_name_fate_file,
)
from name_satellite.prepare import run_settings
from name_satellite.release_met import read_release_met, release_met_soundings
from name_satellite.releases import MODELLED_LEVELS, PARTICLE_RATES

FACES = ("n", "e", "s", "w")
FACE_NAMES = {"n": "N", "e": "E", "s": "S", "w": "W"}
# Boundary heights: m above the model ground, the site files' 20 bins of 1 km (D24).
HEIGHT_BIN_EDGES_M = np.arange(0.0, 20000.0 + 1.0, 1000.0)
# Large arrays are rounded to this many significant bits when written (D23, spec "Storage").
BITROUND_BITS = 10
FATE_LEFT_DOMAIN = 8  # patched NAME fate codes (X20, dev/notes/name-behaviour.md)
FATE_MAX_TRAVEL_TIME = 16
# Half a unit in the last printed place of the fate file's positions (5 decimals; heights 4).
ON_BOUNDARY_TOLERANCE_DEG = 0.5e-5
ON_TOP_TOLERANCE_M = 0.5e-4
SECONDS_PER_HOUR = 3600.0
PA_PER_HPA = 100.0


# ---------------------------------------------------------------------------------------------
# Reading one day's run
# ---------------------------------------------------------------------------------------------


def read_day(day_folder, allow_test=False):
    """Read everything one day's processing needs from its run folder.

    Args:
        day_folder: e.g. ``WORK/GOSAT-BRAZIL_GOSAT-1_v10.0_201301/20130101``.
        allow_test: accept a ``--test`` run (1-hour travel; only to check the code runs).

    Returns:
        dict with ``info`` (day.json), ``site`` (``domains.Site``), ``settings`` (the NAME run
        settings, ``prepare.run_settings``), ``name_version`` (as NAME prints it), ``obs``
        (Dataset keyed by ``sounding_id``), ``members`` (DataFrame: sounding_id -> release_group),
        ``plan``
        (release plan, one row per NAME source), ``fields`` (DataArray (source, lat, lon) in
        ppm s), ``particles`` (fate records, one per particle) and ``met`` (DataFrame, one row per
        release group).

    Raises:
        ValueError if it is a ``--test`` run (unless allowed), the fields are not in ppm s, the
        fate heights are not m agl, or the fate records are not as expected (X20).
        FileNotFoundError if a file is missing.
    """
    day_folder = Path(day_folder)
    info = json.loads((day_folder / "day.json").read_text())
    if info.get("test") and not allow_test:
        raise ValueError(f"{day_folder} is a --test run (1-hour travel): not for footprints")
    site = SITES[info["site"]]
    day = pd.Timestamp(info["day"])
    with xr.open_dataset(day_folder / f"obs_{site.name}_{day:%Y%m%d}.nc") as f:
        obs = f.load()
    members = pd.read_csv(day_folder / "release_groups.csv", dtype=str)
    plan = pd.read_csv(
        day_folder / "release_plan.csv",
        parse_dates=["release_start", "release_end"],
        dtype={"release_group": str},
    )
    if "release_group" not in plan:
        raise ValueError(
            f"{day_folder}/release_plan.csv has no release_group column: prepared before "
            "2026-11-13 (when it was called sounding_id); prepare the day again"
        )
    output = day_folder / "output"
    fields = read_fields_folder(output)
    if set(fields["units"].values) != {"ppm s"}:
        raise ValueError(f"Expected mixing ratio in ppm s, got {set(fields['units'].values)}")
    particles, z_units = read_name_fate_file(output / "Particles_C1.txt")
    if z_units != "m agl":
        raise ValueError(f"Fate heights are in {z_units}; expected m agl (D24)")
    settings = run_settings(site, day_folder, day, info["threads"], info["test"])
    # X20 repeats are one particle time step apart, and a step is never longer than the sync time.
    particles = drop_repeated_fate_records(particles, max_gap_seconds=60.0 * settings.sync_minutes)
    met = read_release_met(output, release_met_soundings(plan))
    # The NAME version is the first line of every NAME output file, e.g. "NAME III (version 8.7.1)".
    fields_file = sorted(output.glob("Fields_grid1*_C1_T1_*.txt*"))[0]
    with open(fields_file) as f:
        name_version = f.readline().strip()
    return {"info": info, "site": site, "obs": obs, "members": members, "plan": plan,
            "fields": fields, "particles": particles, "met": met, "name_version": name_version,
            "settings": settings}  # fmt: skip


# ---------------------------------------------------------------------------------------------
# Particle fates
# ---------------------------------------------------------------------------------------------


def classify_fates(particles, domain):
    """Assign each particle to one fate: "n", "e", "s", "w" (a side face), "top" or "inside".

    NAME records a particle that leaves the domain at its first position outside it (one time
    step beyond the boundary), with fate code 8; one still inside at the maximum travel time has
    code 16. The codes are bit flags: a particle that leaves on its last step has 24 (8 + 16; one
    of 545,000 in a 30-day run, beyond the east edge), and is treated as having left. So, for
    codes 8 and 24:

    * above the domain top (height > ``domain.z_max_m_agl``) -> "top" (also when it is beyond a
      side as well: it cannot go in a side face's height bins, which end at the top);
    * otherwise beyond a side -> that face; beyond two sides at once (a corner) -> the side it is
      further beyond, in degrees;
    * otherwise (code 8 but inside the domain) -> error.

    "Beyond" includes "on the boundary as printed": the fate file gives positions to 5 decimal
    places (about 1 m), and in a 30-day run 31 of 536,000 exits sat exactly on the east edge as
    printed (``ON_BOUNDARY_TOLERANCE_DEG``; heights likewise ``ON_TOP_TOLERANCE_M``). These
    tolerances only place records that are strictly outside nowhere; a record strictly above the
    top or strictly beyond a side is placed by the rules above.

    Args:
        particles: fate records (``read_name_fate_file`` from the patched NAME, one per particle).
        domain: ``name_input.DispersionDomain`` of the run.

    Returns:
        Array of fate names, one per row of ``particles``.

    Raises:
        ValueError for an unknown fate code, a code-8 particle that is inside the domain, or a
        code-16 particle that is outside it.
    """
    codes = particles["fate_code"].to_numpy()
    known = [FATE_LEFT_DOMAIN, FATE_MAX_TRAVEL_TIME, FATE_LEFT_DOMAIN | FATE_MAX_TRAVEL_TIME]
    unknown = ~np.isin(codes, known)
    if unknown.any():
        raise ValueError(f"Unexpected fate codes {sorted(set(codes[unknown]))}")

    lon = particles["lon"].to_numpy()
    lat = particles["lat"].to_numpy()
    z = particles["z"].to_numpy()
    # How far beyond each side (degrees; <= 0 means not beyond it).
    beyond = {
        "n": lat - domain.y_max,
        "e": lon - domain.x_max,
        "s": domain.y_min - lat,
        "w": domain.x_min - lon,
    }
    beyond_table = np.column_stack([beyond[f] for f in FACES])
    furthest = np.array(FACES)[beyond_table.argmax(axis=1)]
    outside_sides = beyond_table.max(axis=1) > 0
    above_top = z > domain.z_max_m_agl
    # For particles NAME says have left: on the boundary (to the printed precision) counts.
    at_or_beyond_sides = beyond_table.max(axis=1) >= -ON_BOUNDARY_TOLERANCE_DEG
    at_or_above_top = z >= domain.z_max_m_agl - ON_TOP_TOLERANCE_M

    left = (codes & FATE_LEFT_DOMAIN) != 0
    # Strictly outside decides first; the tolerances only place records that are outside nowhere.
    fate = np.select(
        [above_top, outside_sides, at_or_above_top], ["top", furthest, "top"], default=furthest
    ).astype(object)
    fate[~left] = "inside"
    unexplained = left & ~at_or_above_top & ~at_or_beyond_sides
    if unexplained.any():
        raise ValueError(f"{unexplained.sum()} particles left the domain but are inside it")
    still_inside = ~left
    if (still_inside & (outside_sides | above_top)).any():
        raise ValueError("Particles that reached the maximum travel time are outside the domain")
    return fate


def _nearest_cell(values, centres):
    """Index of the nearest of the evenly spaced ``centres``, clipped to the grid (D13)."""
    step = (centres[-1] - centres[0]) / (len(centres) - 1)
    index = np.rint((values - centres[0]) / step).astype(int)
    return np.clip(index, 0, len(centres) - 1)


def release_fates(particles, fate, source_numbers, lat, lon):
    """Particle-fate histograms of every release, as fractions of the particles it released.

    Args:
        particles: fate records, one per particle (``source_number``, ``particle_index``,
            ``age_hours``, ``lon``, ``lat``, ``z`` in m agl).
        fate: from ``classify_fates``.
        source_numbers: NAME source numbers of the releases, in the order wanted.
        lat, lon: output-grid cell centres (degrees).

    Returns:
        dict of arrays, first dimension the release (in ``source_numbers`` order):
          * ``n_released`` (release): particles released, counted from the fate records;
          * ``n``, ``s`` (release, height bin, lon) and ``e``, ``w`` (release, height bin, lat):
            fraction of the release's particles leaving through each face, by height bin above
            ground (``HEIGHT_BIN_EDGES_M``) and output-grid cell along the edge (the nearest,
            clipped to the grid);
          * ``top`` (release, lat, lon): fraction leaving through the top, by nearest grid cell;
          * ``inside`` (release): fraction still inside at the maximum travel time;
          * ``age_<fate>`` (same shape as the histogram): sum over those particles of travel time
            (hours) divided by the particles released, so ``age / fraction`` is the mean age.

    Raises:
        ValueError if a release has no records, its particle indices are not 1..N without gaps,
        or a side exit is below the ground or above the top bin.
    """
    source = particles["source_number"].to_numpy()
    wanted = np.asarray(source_numbers)
    position = {number: i for i, number in enumerate(wanted)}
    if not set(np.unique(source)) <= set(wanted):
        raise ValueError("Fate records from sources that are not in the release plan")
    release = np.array([position[number] for number in source])
    n_releases = len(wanted)

    counts = particles.groupby("source_number")["particle_index"].agg(["count", "max", "min"])
    missing = sorted(set(wanted) - set(counts.index))
    if missing:
        raise ValueError(f"No fate records for sources {missing[:5]}")
    counts = counts.loc[wanted]
    if ((counts["min"] != 1) | (counts["max"] != counts["count"])).any():
        raise ValueError("Particle indices are not 1..N for every source (X20)")
    n_released = counts["count"].to_numpy().astype(float)

    age = particles["age_hours"].to_numpy()
    z = particles["z"].to_numpy()
    n_heights = len(HEIGHT_BIN_EDGES_M) - 1
    result = {"n_released": n_released}
    for face in FACES:
        on_face = fate == face
        along_centres = lon if face in ("n", "s") else lat
        along_values = (particles["lon"] if face in ("n", "s") else particles["lat"]).to_numpy()
        heights = z[on_face]
        if len(heights) and (heights.min() < 0 or heights.max() > HEIGHT_BIN_EDGES_M[-1]):
            raise ValueError(f"Side exits outside the height bins on face {face}")
        # Bins [0, 1000), [1000, 2000), ...; a particle exactly at the top goes in the last bin.
        height_bin = np.minimum(
            np.searchsorted(HEIGHT_BIN_EDGES_M, heights, side="right") - 1, n_heights - 1
        )
        along = _nearest_cell(along_values[on_face], along_centres)
        shape = (n_releases, n_heights, len(along_centres))
        result[face] = np.zeros(shape)
        result[f"age_{face}"] = np.zeros(shape)
        index = (release[on_face], height_bin, along)
        np.add.at(result[face], index, 1.0)
        np.add.at(result[f"age_{face}"], index, age[on_face])

    on_top = fate == "top"
    shape = (n_releases, len(lat), len(lon))
    result["top"] = np.zeros(shape)
    result["age_top"] = np.zeros(shape)
    index = (
        release[on_top],
        _nearest_cell(particles["lat"].to_numpy()[on_top], lat),
        _nearest_cell(particles["lon"].to_numpy()[on_top], lon),
    )
    np.add.at(result["top"], index, 1.0)
    np.add.at(result["age_top"], index, age[on_top])

    result["inside"] = np.bincount(release[fate == "inside"], minlength=n_releases).astype(float)

    # Fractions of the particles released (D13): every release's fates now sum to 1.
    for name in [*FACES, "top", "inside", *(f"age_{f}" for f in (*FACES, "top"))]:
        result[name] /= n_released.reshape((-1,) + (1,) * (result[name].ndim - 1))
    return result


# ---------------------------------------------------------------------------------------------
# One day
# ---------------------------------------------------------------------------------------------


def _mean_age(age_sum, fraction):
    """Mean age (hours) from summed ages and fractions; 0 where no particles (as the site files)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(fraction > 0, age_sum / fraction, 0.0)


def process_day(day_folder, allow_test=False):
    """Column footprints, boundary histograms and met for every sounding of one day.

    Args:
        day_folder: one day's run folder (see the module docstring).
        allow_test: accept a ``--test`` run (``read_day``).

    Returns:
        xarray.Dataset in the layout of ``docs/footprint-file-spec.md``, values in float64
        (rounding happens only when writing). Soundings in obs-file order; ``footprint`` equals
        ``sounding_id`` (``footprint_index`` set when the month is assembled).

    Raises:
        ValueError if the inputs disagree (sources, groups, levels) or the particle accounting
        does not add up.
    """
    day = read_day(day_folder, allow_test)
    obs, plan, members, fields = day["obs"], day["plan"], day["members"], day["fields"]
    product = day["info"]["product"]
    n_levels = MODELLED_LEVELS[product]
    domain = day["site"].domain

    # --- Match sources, groups and soundings by their IDs, never by order. -------------------
    plan = plan.set_index("source_name", drop=False)
    if sorted(fields["source"].values) != sorted(plan.index):
        raise ValueError("Fields sources and release plan sources differ")
    plan["source_number"] = [int(name.removeprefix("SourceID_")) for name in plan.index]
    groups = plan.groupby("release_group")
    for group_id, releases in groups:
        if sorted(releases["level_number"]) != list(range(1, n_levels + 1)):
            raise ValueError(f"Group {group_id} does not have release levels 1..{n_levels}")
    sounding_ids = [str(s) for s in obs["sounding_id"].values]
    group_of = members.set_index("sounding_id")["release_group"]
    if sorted(group_of.index) != sorted(sounding_ids):
        raise ValueError("release_groups.csv and the observation file have different soundings")
    if not set(group_of.values) == set(groups.groups):
        raise ValueError("release_groups.csv and the release plan have different groups")

    # Release-plan order: group by group (first appearance), level 1..n within each.
    group_order = list(dict.fromkeys(plan["release_group"]))
    ordered = plan.sort_values(["level_number"]).sort_values(
        "release_group", key=lambda s: s.map(group_order.index), kind="stable"
    )
    source_names = ordered.index.to_numpy().reshape(len(group_order), n_levels)

    # --- 6a: per-release footprints, (group, level, lat, lon). ------------------------------
    lat = fields["lat"].values
    lon = fields["lon"].values
    durations = (plan["release_start"] - plan["release_end"]).dt.total_seconds().unique()
    if len(durations) != 1:
        raise ValueError(f"Releases of different durations: {durations}")
    area = grid_cell_area_m2(lat, lon)
    per_release = footprint_from_mixing_ratio(
        fields.sel(source=source_names.ravel()).values, area, float(durations[0])
    ).reshape(len(group_order), n_levels, len(lat), len(lon))

    # --- 6b: particle fates per release. -------------------------------------------------
    particles = day["particles"]
    fate = classify_fates(particles, domain.dispersion)
    numbers = plan.loc[source_names.ravel(), "source_number"].to_numpy()
    fates = release_fates(particles, fate, numbers, lat, lon)
    level_shape = (len(group_order), n_levels)
    fates = {k: v.reshape(level_shape + v.shape[1:]) for k, v in fates.items()}
    total = sum(fates[f].sum(axis=tuple(range(2, fates[f].ndim))) for f in (*FACES, "top"))
    total = total + fates["inside"]
    if not np.allclose(total, 1.0, rtol=0, atol=1e-12):
        raise ValueError("Particle fates of a release do not sum to 1")

    # --- 6c: columns, one per sounding. -------------------------------------------------
    pressure_weights = obs["pressure_weights"].values[:, :n_levels].astype(float)
    kernel = obs["xch4_averaging_kernel"].values[:, :n_levels].astype(float)
    weights = pressure_weights * kernel  # (sounding, level)
    if not np.isfinite(weights).all():
        raise ValueError("Missing pressure weights or averaging kernels in the modelled levels")
    group_index = np.array([group_order.index(group_of[s]) for s in sounding_ids])

    def column(per_level):
        """Σ_j w_sj × per_level[group of s, j] for every sounding s (one sounding at a time, so
        a busy day's per-level arrays are not copied once per sounding)."""
        sums = [np.tensordot(weights[s], per_level[group_index[s]], axes=1)
                for s in range(len(sounding_ids))]  # fmt: skip
        return np.stack(sums)

    data = {"srr": column(per_release)}
    for face in FACES:
        data[f"particle_locations_{face}"] = column(fates[face])
        data[f"mean_age_particles_{face}"] = _mean_age(
            column(fates[f"age_{face}"]), data[f"particle_locations_{face}"]
        )
    data["particle_locations_top"] = column(fates["top"])
    data["mean_age_particles_top"] = _mean_age(
        column(fates["age_top"]), data["particle_locations_top"]
    )
    data["particle_fraction_inside"] = column(fates["inside"])

    # Every particle counted once: the column fates add up to the column weight.
    column_total = data["particle_fraction_inside"] + data["particle_locations_top"].sum((1, 2))
    for face in FACES:
        column_total = column_total + data[f"particle_locations_{face}"].sum((1, 2))
    if not np.allclose(column_total, weights.sum(axis=1), rtol=1e-12, atol=0):
        raise ValueError("Column particle fates do not add up to the column weight")

    # --- Release layers and met of each sounding's group. ----------------------------------
    release_pressure = plan.loc[source_names.ravel(), "z_pa"].to_numpy().reshape(level_shape)
    release_thickness = plan.loc[source_names.ravel(), "dz_pa"].to_numpy().reshape(level_shape)
    level_1 = ordered[ordered["level_number"] == 1].set_index("release_group").loc[group_order]
    met = day["met"].set_index("release_group").loc[group_order]
    g = group_index
    s_dim, fp_dim = "sounding_id", "footprint"
    rates = PARTICLE_RATES[product]
    ds = xr.Dataset(
        {
            "srr": ((fp_dim, "latitude", "longitude"), data["srr"]),
            **{
                f"particle_locations_{f}": (
                    (fp_dim, "height", "longitude" if f in ("n", "s") else "latitude"),
                    data[f"particle_locations_{f}"],
                )
                for f in FACES
            },
            **{
                f"mean_age_particles_{f}": (
                    (fp_dim, "height", "longitude" if f in ("n", "s") else "latitude"),
                    data[f"mean_age_particles_{f}"],
                )
                for f in FACES
            },
            "particle_locations_top": ((fp_dim, "latitude", "longitude"),
                                       data["particle_locations_top"]),
            "mean_age_particles_top": ((fp_dim, "latitude", "longitude"),
                                       data["mean_age_particles_top"]),
            "particle_fraction_inside": (fp_dim, data["particle_fraction_inside"]),
            "footprint_sounding_id": (fp_dim, np.array(sounding_ids)),
            "time": (s_dim, pd.to_datetime(obs["time"].values)),
            "release_group": (s_dim, np.array([group_of[s] for s in sounding_ids])),
            "air_temperature": (s_dim, met["temperature_k"].to_numpy()[g]),
            "air_pressure": (s_dim, met["pressure_pa"].to_numpy()[g] / PA_PER_HPA),
            "wind_speed": (s_dim, met["wind_speed_m_s"].to_numpy()[g]),
            "wind_from_direction": (s_dim, met["wind_direction_deg"].to_numpy()[g]),
            "atmosphere_boundary_layer_thickness": (
                s_dim, met["boundary_layer_depth_m"].to_numpy()[g]),
            "surface_pressure": (s_dim, met["surface_pressure_pa"].to_numpy()[g] / PA_PER_HPA),
            "release_lat": (s_dim, level_1["lat"].to_numpy()[g]),
            "release_lon": (s_dim, level_1["lon"].to_numpy()[g]),
            "release_dlat": (s_dim, level_1["dlat"].to_numpy()[g]),
            "release_dlon": (s_dim, level_1["dlon"].to_numpy()[g]),
            "column_weight": ((s_dim, "release_level"), weights),
            "release_pressure": ((s_dim, "release_level"), release_pressure[g] / PA_PER_HPA),
            "release_thickness": ((s_dim, "release_level"), release_thickness[g] / PA_PER_HPA),
            "particles_released": ((s_dim, "release_level"), fates["n_released"][g]),
        },
        coords={
            s_dim: np.array(sounding_ids),
            fp_dim: np.arange(len(sounding_ids)),
            "latitude": lat,
            "longitude": lon,
            "height": (HEIGHT_BIN_EDGES_M[:-1] + HEIGHT_BIN_EDGES_M[1:]) / 2.0,
            "release_level": np.arange(1, n_levels + 1),
        },
        attrs={
            "product": product,
            "site": day["site"].name,
            "domain": domain.name,
            "l2_files": obs.attrs.get("source_file", ""),
            "max_level": n_levels,
            "particles_per_second": ", ".join(str(r) for r in rates),
            "model_version": f"{day['name_version']}, ACRG fate-record patch (X20)",
            "met": day["settings"].met_version.name,
            "max_travel_hours": day["settings"].max_age_hours,
            "integration_window": day["settings"].integration,
            "random_seed": day["settings"].random_seed,
        },
    )  # fmt: skip
    return ds


# ---------------------------------------------------------------------------------------------
# A month, and the file
# ---------------------------------------------------------------------------------------------


def process_month(month_folder, days=None):
    """Process every day folder of a month and put the days together (``assemble``).

    Args:
        month_folder: ``WORK/<site>_<product>_<YYYYMM>`` (``prepare.prepare_month``).
        days: optional list of "YYYYMMDD" folders to use (default: every day folder present).

    Returns:
        The month's Dataset (``assemble``).
    """
    month_folder = Path(month_folder)
    if days is None:
        days = sorted(p.name for p in month_folder.iterdir() if p.is_dir() and p.name.isdigit())
    return assemble([process_day(month_folder / d) for d in days])


def assemble(parts):
    """Put days together, sort the soundings and set the footprint index.

    Args:
        parts: ``process_day`` results.

    Returns:
        Dataset of all soundings sorted by time, then sounding ID. Footprints are in the same
        order, so ``footprint_index`` is 0, 1, 2, ... (one exact footprint per sounding, D23);
        ``footprint_sounding_id`` names each footprint's sounding.

    Raises:
        ValueError if a sounding appears twice, or the days differ in product, site, domain,
        modelled levels, particle rates or grid.
    """
    for part in parts[1:]:
        for key in (
            "product",
            "site",
            "domain",
            "max_level",
            "particles_per_second",
            "model_version",
            "max_travel_hours",
            "integration_window",
            "random_seed",
        ):
            if part.attrs[key] != parts[0].attrs[key]:
                raise ValueError(f"Days differ in {key}")
    footprint_vars = [v for v in parts[0].data_vars if "footprint" in parts[0][v].dims]
    sounding_vars = [v for v in parts[0].data_vars if v not in footprint_vars]
    footprints = xr.concat(
        [p[footprint_vars].drop_vars("footprint") for p in parts], dim="footprint", join="exact"
    )
    soundings = xr.concat([p[sounding_vars] for p in parts], dim="sounding_id", join="exact")

    ids = soundings["sounding_id"].values.astype(str)
    if len(set(ids)) != len(ids):
        raise ValueError("A sounding appears more than once")
    order = np.lexsort((ids, soundings["time"].values))
    soundings = soundings.isel(sounding_id=order)
    footprints = footprints.isel(footprint=order)
    if not np.array_equal(footprints["footprint_sounding_id"].values, soundings["sounding_id"]):
        raise ValueError("Footprints and soundings are out of step")

    ds = xr.merge([soundings, footprints], join="exact", combine_attrs="override")
    ds = ds.assign_coords(footprint=np.arange(ds.sizes["footprint"]))
    ds["footprint_index"] = ("sounding_id", np.arange(ds.sizes["sounding_id"], dtype="int32"))
    ds.attrs["l2_files"] = ", ".join(p.attrs["l2_files"] for p in parts)
    ds.attrs["met"] = ", ".join(dict.fromkeys(p.attrs["met"] for p in parts))
    return ds


# ---------------------------------------------------------------------------------------------
# Writing the file (docs/footprint-file-spec.md)
# ---------------------------------------------------------------------------------------------

# Deviations from legacy that new-run files carry (IDs from docs/deviations.md).
DEVIATIONS_APPLIED = "X1 X4 X5 X6 X7 X9 X10 X11 X14 X17 X18 X20 X21 X22 X24 X25"
DECISIONS_APPLIED = "D13 D17 D18 D19 D20 D21 D22 D23 D24"

LOSS_ATTRS = {
    "loss_lifetime_hrs": -9.0,
    "loss_lifetime_comment": "lifetime in hours; -9 corresponds to inert",
}
SRR_UNITS = "(mol/mol)/(mol/m2/s)"

# name -> (units, long_name); names and wording follow the site footprint files where they exist.
VARIABLE_ATTRS = {
    "srr": (SRR_UNITS, "source_receptor_relationship (column)"),
    "particle_locations_top": ("1", "Column-weighted fraction of particles leaving domain (top)"),
    "mean_age_particles_top": ("hrs", "Mean age of particles leaving domain (top)"),
    "particle_fraction_inside": (
        "1",
        "Column-weighted fraction of particles still inside the domain at the maximum travel time",
    ),
    "footprint_sounding_id": ("", "sounding whose column weights built this footprint"),
    "footprint_index": (
        "1",
        "index of each sounding's footprint along the footprint dimension (0-based)",
    ),
    "release_group": ("", "sounding ID of the earliest sounding sharing these NAME releases"),
    "air_temperature": ("K", "air temperature at release (lowest release layer)"),
    "air_pressure": ("hPa", "air pressure at release (lowest release layer)"),
    "wind_speed": ("m s-1", "wind speed at release (lowest release layer)"),
    "wind_from_direction": ("degree", "wind direction at release (lowest release layer)"),
    "atmosphere_boundary_layer_thickness": ("m", "atmospheric boundary layer thickness at release"),
    "surface_pressure": ("hPa", "model surface pressure at the release point"),
    "release_lat": ("degrees_north", "Release latitude (centre of the release box)"),
    "release_lon": ("degrees_east", "Release longitude (centre of the release box)"),
    "release_dlat": ("degrees", "Release box full height"),
    "release_dlon": ("degrees", "Release box full width"),
    "column_weight": ("1", "pressure weight x averaging kernel of each modelled level"),
    "release_pressure": ("hPa", "pressure at the centre of each release layer"),
    "release_thickness": ("hPa", "pressure thickness of each release layer"),
    "particles_released": ("1", "NAME particles released at each level (from the fate records)"),
}
for _f in FACES:
    VARIABLE_ATTRS[f"particle_locations_{_f}"] = (
        "1",
        f"Column-weighted fraction of particles leaving domain ({FACE_NAMES[_f]} side)",
    )
    VARIABLE_ATTRS[f"mean_age_particles_{_f}"] = (
        "hrs",
        f"Mean age of particles leaving domain ({FACE_NAMES[_f]} side)",
    )
ROUNDED = [
    "srr",
    "particle_locations_top",
    "mean_age_particles_top",
    *(f"particle_locations_{f}" for f in FACES),
    *(f"mean_age_particles_{f}" for f in FACES),
]


def software_commit():
    """Git commit of this code (with "-dirty" if it has uncommitted changes), or "unknown"."""
    import subprocess

    try:
        return subprocess.run(
            ["git", "describe", "--always", "--dirty"],
            capture_output=True,
            text=True,
            check=True,
            cwd=Path(__file__).parent,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def write_footprint_file(ds, path):
    """Write a ``process_month``/``assemble`` Dataset as a footprint file (spec, D23).

    Large arrays (``ROUNDED``) are float32 rounded to ``BITROUND_BITS`` significant bits
    (netCDF BitRound: relative precision 2**-(bits+1) = 0.05%, unbiased, small values kept) and
    zlib-compressed, one footprint per chunk. Everything else is stored exactly. The file must
    not exist already.

    Args:
        ds: Dataset from ``assemble``.
        path: output file, e.g. ``GOSAT1v10-BRAZIL-column_NAME_UMG_SOUTHAMERICA_inert_201301.nc``.

    Raises:
        FileExistsError if ``path`` exists.
    """
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    ds = ds.copy()
    for name, (units, long_name) in VARIABLE_ATTRS.items():
        if name in ds:
            ds[name].attrs.update({"units": units, "long_name": long_name})
    ds["srr"].attrs.update(LOSS_ATTRS)
    ds["latitude"].attrs.update({"units": "degrees_north", "long_name": "latitude"})
    ds["longitude"].attrs.update({"units": "degrees_east", "long_name": "longitude"})
    ds["height"].attrs.update(
        {
            "units": "m",
            "positive": "up",
            "bounds": "height_bounds",
            "long_name": "height above model ground at layer midpoints (D24)",
        }
    )
    ds["height_bounds"] = (
        ("height", "bound"),
        np.column_stack([HEIGHT_BIN_EDGES_M[:-1], HEIGHT_BIN_EDGES_M[1:]]),
    )
    ds["release_level"].attrs.update({"long_name": "release level (1 = surface layer)"})
    ds["time"].attrs.update(
        {
            "long_name": "time",
            "standard_name": "time",
            "label": "left",
            "comment": "measurement time from the L2 file: whole minutes, rounded down (X22)",
        }
    )
    month_start = pd.Timestamp(ds["time"].values.min()).to_period("M").start_time
    ds.attrs.update(
        {
            "species": "inert",
            "model": "NAME",
            "met_model": "UMG",
            "lpdm_native_output_unit": "ppm s",
            "author": "University of Bristol",
            "publication_acknowledgement": "Please acknowledge the Met Office/University of Bristol in any publications that use these files.",  # noqa: E501
            "deviations_applied": DEVIATIONS_APPLIED,
            "decisions_applied": DECISIONS_APPLIED,
            "software": "name-satellite",
            "software_commit": software_commit(),
            "created": pd.Timestamp.now(tz="UTC").isoformat(timespec="seconds"),
            "Conventions": "CF-1.8",
        }
    )

    encoding = {
        "time": {
            "units": f"seconds since {month_start:%Y-%m-%d %H:%M:%S}",
            "calendar": "gregorian",
            "dtype": "float64",
        },
        "footprint_index": {"dtype": "int32"},
    }
    for name in ROUNDED:
        chunks = (1, *ds[name].shape[1:])
        encoding[name] = {
            "dtype": "float32",
            "zlib": True,
            "complevel": 4,
            "shuffle": True,
            "chunksizes": chunks,
            "quantize_mode": "BitRound",
            "significant_digits": BITROUND_BITS,
        }
    for name in ds.data_vars:
        if name not in encoding and ds[name].dtype.kind == "f":
            encoding[name] = {"dtype": "float64", "zlib": True}
    ds.to_netcdf(path, engine="netcdf4", format="NETCDF4", encoding=encoding)


def footprint_file_name(ds):
    """File name in the site files' pattern (``{site}_{model}_{met}_{domain}_{species}_{YYYYMM}``).

    The site part (product, version, region, "column") is still to be agreed with the group
    (spec open question), e.g. ``GOSAT1v10-BRAZIL-column_NAME_UMG_SOUTHAMERICA_inert_201301.nc``.
    """
    product = ds.attrs["product"].replace("GOSAT-", "GOSAT").replace(" v", "v").split(".")[0]
    region = ds.attrs["site"].split("-", 1)[-1]
    month = pd.Timestamp(ds["time"].values.min())
    return f"{product}-{region}-column_NAME_UMG_{ds.attrs['domain']}_inert_{month:%Y%m}.nc"


if __name__ == "__main__":
    import sys

    usage = "uv run python -m name_satellite.processing MONTH_FOLDER OUT_FOLDER [--days=YYYYMMDD]"
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 2:
        sys.exit(usage)
    chosen = [a.removeprefix("--days=").split(",") for a in sys.argv[1:] if a.startswith("--days=")]
    month = process_month(args[0], days=chosen[0] if chosen else None)
    Path(args[1]).mkdir(parents=True, exist_ok=True)
    out = Path(args[1]) / footprint_file_name(month)
    write_footprint_file(month, out)
    print(f"{out}: {month.sizes['sounding_id']} soundings")
