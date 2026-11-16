"""Turn NAME output into the files GLIDE asks for (spec §8): footprints and end records.

Where this sits
---------------
    prepare.py  ->  LOTUS2 (NAME)  ->  **deliver.py: one footprint file and one end-record file
    per receptor, block and run**

Inputs: one or more run folders (``prepare.py``) that together hold the releases of one
receptor, block and run variant (a block may be split into several runs, e.g. for memory). From
each folder: ``plan.csv``, ``run.json``, the fields files and the particle-fate records
(``Particles_C1.txt``, ACRG-patched NAME) in ``output/``.

Outputs (spec §8; names ``{receptor}_NAME_UMG-Mk10_{domain}_{YYYYMM}_{run}.nc`` and ``..._end.nc``):

* **Footprint file:** ``srr(time, latitude, longitude)``, float32, (mol/mol)/(mol/m2/s); ``time``
  is the release start; ``release_lon``, ``release_lat`` and ``release_height`` (m above the
  model's ground) per release, or for the column ``release_pressure_bottom`` and
  ``release_pressure_top`` (hPa); ``n_particles`` (counted from the fate records);
  ``fraction_left_sides`` and ``fraction_left_top`` (of ``n_particles``).
* **End-record file:** exactly one row per particle (dimension ``particle``): its release
  (``time``), ``particle_index``, the position where NAME stopped following it (``lon``, ``lat``,
  ``height`` m above the model's ground), ``age`` (hours) and ``fate``: "n", "e", "s", "w" (left
  through that side), "top" (left through the top) or "max_age" (reached 240 h inside the
  domain). ``reached_max_age`` is also true for a particle that left on the step it reached 240 h
  (NAME gives both reasons; §8 asks for all).

The footprint (spec §6) is NAME's time-integrated mixing ratio in the 0-40 m layer, converted as
in the group's footprints (``footprints.footprint_from_mixing_ratio``): field (ppm s) × 1e-6 ×
cell area (m², sphere of radius 6,367,500 m) / release duration (3,600 s). NAME divides each
source's mass among the particles it actually released, so the footprint is already normalised by
the true number of particles (§4, §6).

Releases are matched to NAME's output by source name (``SourceID_<n>``, the field request
``Req_H1Z1_<n>``) and, for the fate records, by NAME's source number n, which is the source's
position in the input file: ``plan.csv`` lists the sources in that order, which is checked.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from name_satellite.footprints import footprint_from_mixing_ratio, grid_cell_area_m2
from name_satellite.glide_runs import spec
from name_satellite.name_output import (
    drop_repeated_fate_records,
    read_fields_folder,
    read_name_fate_file,
)
from name_satellite.processing import FATE_MAX_TRAVEL_TIME, classify_fates

MODEL = "NAME"
METEOROLOGY = "UMG-Mk10"
RELEASE_SECONDS = spec.RELEASE_HOURS * 3600.0


def file_names(receptor, block, run):
    """(footprint file, end-record file) names for a receptor code, ``spec.Block`` and run name."""
    stem = f"{receptor}_{MODEL}_{METEOROLOGY}_{block.domain}_{block.month}_{run}"
    return f"{stem}.nc", f"{stem}_end.nc"


def read_plan(run_folder):
    """The run's release plan, checked to be in NAME's source order (SourceID_1, 2, ...)."""
    plan = pd.read_csv(
        Path(run_folder) / "plan.csv", parse_dates=["release_time", "release_start", "release_end"]
    )
    expected = [f"SourceID_{i}" for i in range(1, len(plan) + 1)]
    if plan["source_name"].tolist() != expected:
        raise ValueError(f"{run_folder}/plan.csv is not in NAME source order")
    return plan


def end_records(run_folder, plan, domain):
    """One row per particle: where and why NAME stopped following it (see the module docstring).

    Raises:
        ValueError if a release's particles are not numbered 1..N without gaps (a particle
        without a fate record), or for fate records ``processing.classify_fates`` rejects.
    """
    fates, z_units = read_name_fate_file(Path(run_folder) / "output" / "Particles_C1.txt")
    if z_units != "m agl":
        raise ValueError(f"Fate heights must be m agl, not {z_units}")
    fates = drop_repeated_fate_records(fates, max_gap_seconds=spec.SYNC_MINUTES * 60.0)
    fate = classify_fates(fates, domain)
    records = pd.DataFrame(
        {
            "source_number": fates["source_number"].to_numpy(),
            "particle_index": fates["particle_index"].to_numpy(),
            "lon": fates["lon"].to_numpy(),
            "lat": fates["lat"].to_numpy(),
            "height": fates["z"].to_numpy(),
            "age": fates["age_hours"].to_numpy(),
            "fate": np.where(fate == "inside", "max_age", fate),
            "reached_max_age": (fates["fate_code"].to_numpy() & FATE_MAX_TRAVEL_TIME) != 0,
        }
    )
    if not records["source_number"].between(1, len(plan)).all():
        raise ValueError(f"Fate records for sources not in {run_folder}/plan.csv")
    counts = records.groupby("source_number")["particle_index"].agg(["size", "max", "min"])
    counts = counts.reindex(range(1, len(plan) + 1))
    if counts.isna().any().any() or (counts["size"] != counts["max"]).any():
        raise ValueError(f"{run_folder}: some particles have no fate record")
    times = plan["release_time"].to_numpy()
    records.insert(0, "time", times[records["source_number"].to_numpy() - 1])
    records = records.sort_values(["time", "particle_index"], kind="stable")
    return records.drop(columns="source_number").reset_index(drop=True)


def footprints(run_folder, plan):
    """Footprints (mol/mol)/(mol/m2/s), DataArray (time, latitude, longitude), one per release."""
    fields = read_fields_folder(Path(run_folder) / "output")
    if fields["source"].values.tolist() != plan["source_name"].tolist():
        raise ValueError(f"{run_folder}: fields and plan have different sources")
    if set(fields["units"].values.tolist()) != {"ppm s"}:
        raise ValueError(f"{run_folder}: expected fields in ppm s")
    area = grid_cell_area_m2(fields["lat"].values, fields["lon"].values)
    srr = footprint_from_mixing_ratio(fields.values, area, RELEASE_SECONDS)
    return xr.DataArray(
        srr.astype(np.float32),
        dims=("time", "latitude", "longitude"),
        coords={
            "time": plan["release_time"].to_numpy(),
            "latitude": fields["lat"].values,
            "longitude": fields["lon"].values,
        },
    )


def deliver(run_folders, receptor_code, block_name, run, out_folder):
    """Write the footprint and end-record files of one receptor, block and run.

    Args:
        run_folders: the ``prepare.py`` run folders holding all releases of this file.
        receptor_code, block_name: keys of ``spec.RECEPTORS`` and ``spec.BLOCKS``.
        run: run name in the file name: "main", "seed1" ... "seed5", "nomeso", "group".
        out_folder: where to write (must exist).

    Returns:
        (footprint path, end-record path).

    Raises:
        ValueError if a folder is for another receptor, a release time appears twice, or a
        check in ``footprints`` / ``end_records`` fails.
    """
    receptor = spec.RECEPTORS[receptor_code]
    block = spec.BLOCKS[block_name]
    domain = spec.DOMAINS[block.domain]
    plans, srrs, ends = [], [], []
    for folder in run_folders:
        plan = read_plan(folder)
        if set(plan["receptor"]) != {receptor_code}:
            raise ValueError(f"{folder} is not a run for {receptor_code}")
        ends.append(end_records(folder, plan, domain.dispersion))
        srrs.append(footprints(folder, plan))
        plans.append(plan)
    plan = pd.concat(plans, ignore_index=True)
    if plan["release_time"].duplicated().any():
        raise ValueError("A release appears in more than one run folder")
    order = np.argsort(plan["release_time"].to_numpy())
    plan = plan.iloc[order].reset_index(drop=True)
    srr = xr.concat(srrs, dim="time").sortby("time")
    end = pd.concat(ends, ignore_index=True).sort_values(["time", "particle_index"])

    # Per-release counts from the end records.
    by_release = end.groupby("time")
    n_particles = by_release.size().reindex(plan["release_time"]).to_numpy()
    side = end["fate"].isin(["n", "e", "s", "w"])
    top = end["fate"] == "top"
    left_sides = side.groupby(end["time"]).sum().reindex(plan["release_time"]).to_numpy()
    left_top = top.groupby(end["time"]).sum().reindex(plan["release_time"]).to_numpy()

    time = plan["release_time"].to_numpy()
    ds = xr.Dataset({"srr": srr})
    ds["release_lon"] = ("time", plan["lon"].to_numpy())
    ds["release_lat"] = ("time", plan["lat"].to_numpy())
    if receptor.is_column:
        ds["release_pressure_bottom"] = ("time", plan["pressure_bottom_pa"].to_numpy() / 100.0)
        ds["release_pressure_top"] = ("time", plan["pressure_top_pa"].to_numpy() / 100.0)
    else:
        ds["release_height"] = ("time", plan["z"].to_numpy())
    ds["n_particles"] = ("time", n_particles.astype(np.int32))
    ds["fraction_left_sides"] = ("time", left_sides / n_particles)
    ds["fraction_left_top"] = ("time", left_top / n_particles)
    ds = ds.assign_coords(time=time)
    _add_attributes(ds, receptor, block, run)

    end_ds = xr.Dataset.from_dataframe(end.reset_index(drop=True).rename_axis("particle"))
    end_ds["fate"] = end_ds["fate"].astype(str)
    end_ds["lon"].attrs = {"units": "degrees_east"}
    end_ds["lat"].attrs = {"units": "degrees_north"}
    end_ds["height"].attrs = {"units": "m", "long_name": "height above the model's ground"}
    end_ds["age"].attrs = {"units": "hours", "long_name": "travel time when NAME stopped"}
    end_ds["time"].attrs = {"long_name": "start of the particle's release (UTC)"}
    end_ds["fate"].attrs = {
        "long_name": "why NAME stopped following the particle",
        "values": "n, e, s, w: left through that side; top: left through the top (20,000 m "
        "agl); max_age: reached the maximum travel time inside the domain",
    }
    end_ds.attrs = dict(ds.attrs)

    out_folder = Path(out_folder)
    footprint_name, end_name = file_names(receptor_code, block, run)
    encoding = {"srr": {"zlib": True, "complevel": 4}}
    ds.to_netcdf(out_folder / footprint_name, encoding=encoding)
    end_ds.to_netcdf(out_folder / end_name)
    return out_folder / footprint_name, out_folder / end_name


def _add_attributes(ds, receptor, block, run):
    ds["srr"].attrs = {
        "units": "(mol/mol)/(mol/m2/s)",
        "long_name": "footprint: mole-fraction sensitivity to a surface flux, 0-40 m layer, "
        "integrated over 240 h",
    }
    ds["latitude"].attrs = {"units": "degrees_north", "long_name": "cell centre"}
    ds["longitude"].attrs = {"units": "degrees_east", "long_name": "cell centre"}
    ds["time"].attrs = {"long_name": "start of the one-hour release (UTC)"}
    ds["release_lon"].attrs = {"units": "degrees_east"}
    ds["release_lat"].attrs = {"units": "degrees_north"}
    if "release_height" in ds:
        ds["release_height"].attrs = {"units": "m", "long_name": "above the model's ground"}
    else:
        for name in ("release_pressure_bottom", "release_pressure_top"):
            ds[name].attrs = {"units": "hPa"}
    ds["n_particles"].attrs = {"long_name": "particles released, counted from NAME's output"}
    ds["fraction_left_sides"].attrs = {"long_name": "fraction removed through a side edge"}
    ds["fraction_left_top"].attrs = {"long_name": "fraction removed through the top (20 km agl)"}
    ds.attrs = {
        "title": f"NAME reference footprints for GLIDE: {receptor.code}, {block.name}, {run}",
        "model": "NAME 8.7.1 (ACRG-patched fate records), Met Office",
        "meteorology": "Met Office UM global analyses, UMG_Mk10PT",
        "specification": spec.SPEC_VERSION,
        "unit_conversion": "NAME mixing ratio (ppm s, 1 g/s release of a 1 g/mol tracer) x 1e-6 x "
        "cell area (m2, sphere radius 6,367,500 m) / 3,600 s",
        "receptor": receptor.code,
        "site": receptor.site,
        "domain": block.domain,
        "block": block.name,
        "run": run,
        "licence": "NAME output: restricted (Met Office licence); not for public distribution",
    }


# ---------------------------------------------------------------------------------------------
# Particle positions (single-release runs, spec §9)
# ---------------------------------------------------------------------------------------------
POSITION_HOURS_BEFORE_START = (6, 24, 72, 120)  # §9


def read_positions(path):
    """Read NAME's particle-position file (a particle output with a time grid).

    Returns:
        DataFrame: ``source`` (name), ``particle_index``, ``time`` (UTC), ``age`` (hours), ``lon``,
        ``lat``, ``height`` (m above the model's ground). The file gives travel time in seconds.

    Raises:
        ValueError if the expected columns are missing (the heading row is found by its first
        column, "Puff?").
    """
    with open(path) as f:
        heading = next(i for i, line in enumerate(f) if line.lstrip().startswith("Puff?"))
    raw = pd.read_csv(path, skiprows=heading, skipinitialspace=True)
    raw.columns = [c.strip() for c in raw.columns]
    needed = [
        "PP Index",
        "Source",
        "Time",
        "Travel Time",
        "X (Lat-Long)",
        "Y (Lat-Long)",
        "Z (m agl)",
    ]
    missing = [c for c in needed if c not in raw.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    return pd.DataFrame(
        {
            "source": raw["Source"].str.strip(),
            "particle_index": raw["PP Index"].astype(int),
            "time": pd.to_datetime(raw["Time"].str.strip(), format="%d/%m/%Y %H:%M:%S UTC"),
            "age": raw["Travel Time"].astype(float) / 3600.0,
            "lon": raw["X (Lat-Long)"].astype(float),
            "lat": raw["Y (Lat-Long)"].astype(float),
            "height": raw["Z (m agl)"].astype(float),
        }
    )


def particle_positions(run_folder, end):
    """Every particle's position at 6, 24, 72 and 120 h before its release's start (§9).

    Args:
        run_folder: a single-release run (one source, ``prepare.py --positions-hours 241``).
        end: that run's end records (``end_records``), to list every particle, including those
            that had already left the domain (NAME no longer writes their position).

    Returns:
        DataFrame, one row per particle and time: ``hours_before_start``, ``particle_index``,
        ``in_domain`` (False once the particle has left), and ``lon``, ``lat``, ``height``, ``age``
        (NaN when not in the domain).
    """
    plan = read_plan(run_folder)
    if len(plan) != 1:
        raise ValueError(f"{run_folder}: particle positions are for single-release runs")
    start = plan["release_time"].iloc[0]
    positions = read_positions(Path(run_folder) / "output" / "Positions_C1.txt")
    parts = []
    for hours in POSITION_HOURS_BEFORE_START:
        at = positions[positions["time"] == start - pd.Timedelta(hours=hours)]
        if at.empty:
            raise ValueError(f"{run_folder}: no positions {hours} h before the release start")
        table = (
            end[["particle_index"]]
            .merge(at.drop(columns=["source", "time"]), on="particle_index", how="left")
            .assign(hours_before_start=hours)
        )
        table["in_domain"] = table["lon"].notna()
        parts.append(table)
    return pd.concat(parts, ignore_index=True)


def deliver_single(run_folder, receptor_code, block_name, run, out_folder):
    """A single-release run (§9): footprint and end-record files, plus ``..._positions.nc``.

    ``run`` is the name used in the file names. The spec gives none for the single-release runs;
    we use ``case<n>`` (n = the standard case number of §9) until it does.
    """
    footprint_path, end_path = deliver([run_folder], receptor_code, block_name, run, out_folder)
    plan = read_plan(run_folder)
    end = end_records(run_folder, plan, spec.DOMAINS[spec.BLOCKS[block_name].domain].dispersion)
    table = particle_positions(run_folder, end)
    ds = xr.Dataset.from_dataframe(table.set_index(["hours_before_start", "particle_index"]))
    ds["lon"].attrs = {"units": "degrees_east"}
    ds["lat"].attrs = {"units": "degrees_north"}
    ds["height"].attrs = {"units": "m", "long_name": "height above the model's ground"}
    ds["age"].attrs = {"units": "hours", "long_name": "travel time"}
    ds["hours_before_start"].attrs = {"long_name": "hours before the release's start time"}
    ds.attrs = dict(xr.open_dataset(footprint_path).attrs)
    ds.attrs["release_time"] = str(plan["release_time"].iloc[0])
    positions_path = Path(str(footprint_path).removesuffix(".nc") + "_positions.nc")
    ds.to_netcdf(positions_path)
    return footprint_path, end_path, positions_path


def main(argv=None):
    """Command line: ``deliver OUT_FOLDER RECEPTOR BLOCK RUN RUN_FOLDER [RUN_FOLDER ...]``.

    With ``--single``, one single-release run folder, also writing the positions file.
    """
    import argparse

    parser = argparse.ArgumentParser(description="Write the GLIDE §8 files from NAME runs.")
    parser.add_argument("out_folder")
    parser.add_argument("receptor", choices=sorted(spec.RECEPTORS))
    parser.add_argument("block", choices=sorted(spec.BLOCKS))
    parser.add_argument("run")
    parser.add_argument("run_folders", nargs="+")
    parser.add_argument("--single", action="store_true")
    args = parser.parse_args(argv)
    if args.single:
        if len(args.run_folders) != 1:
            raise SystemExit("--single takes one run folder")
        paths = deliver_single(
            args.run_folders[0], args.receptor, args.block, args.run, args.out_folder
        )
    else:
        paths = deliver(args.run_folders, args.receptor, args.block, args.run, args.out_folder)
    for path in paths:
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
