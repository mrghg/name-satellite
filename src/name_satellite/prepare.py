"""Prepare a month of NAME runs from Level-2 satellite data: workflow steps 1-4 chained together.

Where this sits (see docs/workflow.md)
--------------------------------------
    L2 files (CEDA)  ->  **prepare**: per day, observation file + surface-pressure run + job
    ->  on LOTUS2: met restore -> surface-pressure run -> main job (writes the NAME input, runs
    NAME)  ->  read and process the output (steps 5-6)

Each day has TWO NAME runs, because the release layers are scaled to the model's surface pressure
(decision D19), which only a NAME run can report:

``prepare_day`` (now, on the sci node):
  1. downloads the day's L2 file (``l2_files``; days without a file are skipped and reported);
  2. reads it and keeps good soundings in the site's box (``observations``: quality flag, defined
     pressure levels, pressures decreasing from the surface, lat/lon bounds), binned if the site
     says so;
  3. writes the observation file, keyed by sounding ID (D18), and ``day.json`` (what was asked);
  4. writes the met-only surface-pressure run (``surface_pressure``) and its job;
  5. writes the main job, which first runs ``write_main_input`` and then NAME.

``write_main_input`` (in the main job, once the surface-pressure run has finished):
  6. reads the model surface pressure and checks it against the retrieval's
     (``surface_pressure.csv``);
  7. scales every sounding's levels to its model surface pressure (D19), checks that the
     modelled levels (``releases.MODELLED_LEVELS``) lie inside the domain and carry the column
     (D21), groups the soundings
     that share releases (``release_groups``, D20; ``release_groups.csv`` maps each sounding to
     its group) and plans one set of releases per group (``releases.build_release_plan``,
     ``release_plan.csv``), with each level's particle rate from ``releases.PARTICLE_RATES``
     (D22);
  8. writes the NAME input (with release-point met, X21).

``prepare_month`` also writes one met-restore job for the whole month and ``submit.sh``, which
submits the restore job, then each day's surface-pressure job after it and the day's main job
after that (Slurm ``afterok``: a job starts only if the one before succeeded).

Run folder layout (``WORK/<site>_<product>_<YYYYMM>/``)::

    restore_met.sbatch, submit.sh, summary.csv
    YYYYMMDD/obs_<site>_<YYYYMMDD>.nc, day.json, job.sbatch, output/
    YYYYMMDD/surface_pressure/soundings.csv, input.txt, job.sbatch, output/
    written by the main job: YYYYMMDD/surface_pressure.csv, release_groups.csv, release_plan.csv,
                             input.txt

    uv run python -m name_satellite.prepare "GOSAT-1 v10.0" GOSAT-BRAZIL 2013-01 WORK [--test]
        [--days=2013-01-01,2013-01-02]
    uv run python -m name_satellite.prepare main-input WORK/.../YYYYMMDD   (what the job runs)

``--test``: particles travel 1 hour only, debug queue, to check a set-up quickly (lesson of X20).
``--days``: prepare only these days of the month (e.g. one day for a test).

The main run uses the ACRG-patched NAME 8.7.1 (``met.PATCHED_NAME_DIR``), whose fate records carry
each particle's index and fate code (X20), and asks for met columns (incl. pressure, Pa) in the
fate records (D14). The met-only run uses the stock NAME 8.7.1.
"""

import dataclasses
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from name_satellite.domains import SITES
from name_satellite.l2_files import L2FileMissing, download_l2
from name_satellite.lotus import REPO, write_job_script, write_met_restore_job
from name_satellite.met import PATCHED_NAME_DIR, met_version_for_period
from name_satellite.name_input import RunSettings, run_timing, write_name_input
from name_satellite.observations import (
    bad_pressure_mask,
    latlon_mask,
    legacy_bin,
    level_order_mask,
    quality_mask,
    read_l2,
    select,
    to_obs,
    write_obs_file,
)
from name_satellite.release_groups import group_soundings, release_points
from name_satellite.release_met import MAX_SOUNDINGS, release_met_soundings
from name_satellite.releases import (
    MODELLED_LEVELS,
    PARTICLE_RATES,
    build_release_plan,
    check_modelled_levels,
    scale_to_model_surface,
)
from name_satellite.surface_pressure import (
    read_surface_pressure,
    surface_pressure_input_text,
    surface_pressure_settings,
    surface_pressure_table,
)

SCRATCH = Path("/work/scratch-pw5/mattrigby/name-satellite")
L2_CACHE = SCRATCH / "l2_cache"
MET_ROOT = SCRATCH / "met"
THREADS = 4  # about 80% efficient (benchmark, dev/notes/name-behaviour.md)


def select_soundings(table, site):
    """The good soundings of one day's table in the site's box (binned if the site says so)."""
    good = quality_mask(table) & bad_pressure_mask(table) & level_order_mask(table)
    table = select(table, good)
    if site.bin_size is not None:
        # legacy_bin applies the box itself (soundings outside the bins are dropped).
        return legacy_bin(
            table, site.lat_bounds, site.lon_bounds, site.bin_size, footprint_km=site.footprint_km
        )
    return select(table, latlon_mask(table, site.lat_bounds, site.lon_bounds))


def run_settings(site, run_folder, day, threads, test):
    """NAME settings for one day's run (legacy defaults apart from the deviations adopted)."""
    met_version = met_version_for_period(day - pd.Timedelta(days=31), day + pd.Timedelta(days=1))
    settings = RunSettings(
        run_name=run_folder.name,
        output_folder=str(run_folder / "output"),
        grid=site.domain.grid,
        domain=site.domain.dispersion,
        met_version=met_version,
        met_folder=str(MET_ROOT / met_version.name),
        threads=threads,
        fate_met=True,  # pressure (Pa) of each particle's fate, for the D14 pressure bins
    )
    if test:
        settings = dataclasses.replace(settings, max_age_hours=1)
    return settings


def release_boxes(obs):
    """(dlon, dlat) per sounding for binned soundings, else None (the footprint box is used)."""
    return (obs["dlon"].values, obs["dlat"].values) if "dlat" in obs else None


def grouped_release_plan(obs, site, n_levels, pressure_levels_hpa=None, particle_rates=None):
    """Release groups (D20) and their release plan for one day.

    Args:
        obs: the day's observation Dataset (one entry per sounding).
        site: the ``domains.Site`` (release cells, footprint size, output grid).
        n_levels: number of modelled levels (``releases.MODELLED_LEVELS``, D21).
        pressure_levels_hpa: optional (sounding, level) levels to release at (D19-scaled).
        particle_rates: particles/s per level (``releases.PARTICLE_RATES``, D22); default legacy.

    Returns:
        (plan, members): the release plan (one set of releases per group; its ``sounding_id`` is
        the group's ID) and the sounding -> group table.
    """
    groups = group_soundings(obs, site.release_cell, site.domain.grid, site.release_group_seconds)
    points, members = release_points(obs, groups, site.footprint_km, pressure_levels_hpa)
    boxes = (points["dlon"].values, points["dlat"].values)
    plan = build_release_plan(points, n_levels=n_levels, footprint_diameter_km=site.footprint_km,
                              box_degrees=boxes, particle_rates=particle_rates)  # fmt: skip
    return plan, members


def prepare_day(product, site, day, month_folder, threads=THREADS, test=False):
    """Prepare one day: observation file, surface-pressure run and the main job.

    Returns a summary dict (one row of summary.csv).
    """
    day = pd.Timestamp(day)
    row = {"day": f"{day:%Y-%m-%d}", "status": "", "l2_soundings": 0, "soundings": 0, "sources": 0}
    try:
        l2_path = download_l2(product, day, L2_CACHE)
    except L2FileMissing:
        row["status"] = "no L2 file"
        return row
    table = read_l2(l2_path)
    row["l2_soundings"] = table.sizes["sounding"]
    table = select_soundings(table, site)
    row["soundings"] = table.sizes["sounding"]
    if row["soundings"] == 0:
        row["status"] = "no soundings in the box"
        return row
    if row["soundings"] > MAX_SOUNDINGS:
        # Release-point met allows at most MAX_SOUNDINGS per run (release_met.py).
        raise ValueError(f"{day:%Y-%m-%d}: {row['soundings']} soundings > {MAX_SOUNDINGS}")

    run_folder = month_folder / f"{day:%Y%m%d}"
    (run_folder / "output").mkdir(parents=True)
    write_obs_file(table, run_folder / f"obs_{site.name}_{day:%Y%m%d}.nc")
    info = {"product": product, "site": site.name, "day": f"{day:%Y-%m-%d}", "threads": threads,
            "test": test}  # fmt: skip
    (run_folder / "day.json").write_text(json.dumps(info, indent=2) + "\n")

    # Provisional plans from the retrieval's own levels: times and positions are final (only the
    # pressures change, D19). Per sounding: where the met-only run samples. Grouped (D20): the
    # number of releases and the run timing.
    obs = to_obs(table)
    n_levels = MODELLED_LEVELS[product]
    sounding_plan = build_release_plan(obs, n_levels=n_levels,
                                       footprint_diameter_km=site.footprint_km,
                                       box_degrees=release_boxes(obs))  # fmt: skip
    plan, members = grouped_release_plan(obs, site, n_levels,
                                         particle_rates=PARTICLE_RATES[product])  # fmt: skip
    settings = run_settings(site, run_folder, day, threads, test)
    qos = "debug" if test else "standard"

    # The met-only run: model surface pressure at each sounding.
    folder = run_folder / "surface_pressure"
    (folder / "output").mkdir(parents=True)
    soundings = release_met_soundings(sounding_plan)
    soundings.to_csv(folder / "soundings.csv", index=False)
    sp_settings = surface_pressure_settings(settings, folder / "output")
    (folder / "input.txt").write_text(surface_pressure_input_text(soundings, sp_settings))
    write_job_script(folder / "job.sbatch", folder, folder / "input.txt", threads=1, qos=qos,
                     time="00:30:00", memory="8G",
                     job_name=f"{site.name}_{day:%Y%m%d}_ps")  # fmt: skip

    # The main job: write the NAME input from the surface pressure, then run NAME.
    write_main = (
        f"cd {REPO}\n"
        f"~/.local/bin/uv run --no-sync python -m name_satellite.prepare main-input {run_folder}"
        " || exit 1"
    )
    qos, time = ("debug", "00:30:00") if test else ("high", "12:00:00")
    write_job_script(run_folder / "job.sbatch", run_folder, run_folder / "input.txt",
                     threads=threads, qos=qos, time=time, memory="32G",
                     job_name=f"{site.name}_{day:%Y%m%d}", before_name=write_main,
                     name_dir=PATCHED_NAME_DIR)  # fmt: skip
    timing = run_timing(plan, settings)  # unaffected by the pressures
    row.update(
        status="prepared",
        release_groups=members["release_group"].nunique(),
        sources=len(plan),
        met_version=settings.met_version.name,
        met_start=f"{timing['window_start']:%Y-%m-%dT%H:%M}",
        met_end=f"{plan['release_start'].max():%Y-%m-%dT%H:%M}",
    )
    return row


def write_main_input(run_folder):
    """Write a day's release groups, release plan and NAME input from its surface-pressure run.

    Run by the day's main job after the surface-pressure job (``submit.sh``); can be run by hand.

    Writes ``surface_pressure.csv`` (retrieval and model surface pressure, hPa, and their ratio),
    ``release_groups.csv`` (each sounding's release group, D20), ``release_plan.csv`` (one set of
    releases per group, levels scaled to the model surface, D19) and ``input.txt``.

    Raises:
        ValueError if the surface-pressure output is missing, or implausible
        (``surface_pressure.surface_pressure_table``).
    """
    run_folder = Path(run_folder)
    info = json.loads((run_folder / "day.json").read_text())
    site = SITES[info["site"]]
    day = pd.Timestamp(info["day"])
    with xr.open_dataset(run_folder / f"obs_{site.name}_{day:%Y%m%d}.nc") as f:
        obs = f.load()

    folder = run_folder / "surface_pressure"
    soundings = pd.read_csv(folder / "soundings.csv", parse_dates=["time"],
                            dtype={"sounding_id": str})  # fmt: skip
    model = read_surface_pressure(folder / "output", soundings)
    surface = surface_pressure_table(obs, model)
    surface.to_csv(run_folder / "surface_pressure.csv", index=False)

    levels = obs["pressure_levels"].values
    scaled = np.array([scale_to_model_surface(levels[p], surface["model_hpa"].iloc[p])
                       for p in range(len(surface))])  # fmt: skip
    n_levels = MODELLED_LEVELS[info["product"]]
    check_modelled_levels(obs, scaled, surface["domain_top_hpa"].to_numpy(), n_levels)
    rates = PARTICLE_RATES[info["product"]]
    plan, members = grouped_release_plan(obs, site, n_levels, pressure_levels_hpa=scaled,
                                         particle_rates=rates)  # fmt: skip
    members.to_csv(run_folder / "release_groups.csv", index=False)
    plan.to_csv(run_folder / "release_plan.csv", index=False)
    settings = run_settings(site, run_folder, day, info["threads"], info["test"])
    write_name_input(run_folder / "input.txt", plan, settings)
    print(f"{run_folder}: {len(surface)} soundings in {members['release_group'].nunique()} release "
          f"groups, {len(plan)} releases; model/retrieval surface pressure "
          f"{surface['ratio'].min():.4f}-{surface['ratio'].max():.4f}")  # fmt: skip
    return plan


def prepare_month(product, site_name, month, work, threads=THREADS, test=False, days=None):
    """Prepare every day of ``month`` ("YYYY-MM"), or only ``days`` (dates in the month).

    Returns the summary table.
    """
    site = SITES[site_name]
    month = pd.Period(month, freq="M")
    month_folder = Path(work) / f"{site.name}_{product.replace(' ', '_')}_{month:%Y%m}"
    if month_folder.exists():
        raise SystemExit(f"{month_folder} already exists; choose a new folder")
    month_folder.mkdir(parents=True)

    all_days = pd.date_range(month.start_time, month.end_time.normalize(), freq="D")
    if days is not None:
        chosen = pd.DatetimeIndex(pd.to_datetime(list(days))).normalize()
        if not chosen.isin(all_days).all():
            raise ValueError(f"Days {list(days)} are not all in {month}")
        all_days = chosen
    rows = [prepare_day(product, site, day, month_folder, threads, test) for day in all_days]
    summary = pd.DataFrame(rows)
    summary.to_csv(month_folder / "summary.csv", index=False)
    prepared = summary[summary["status"] == "prepared"]
    if prepared.empty:
        print(summary.to_string(index=False))
        raise SystemExit("No day had soundings to run")

    versions = prepared["met_version"].unique()
    if len(versions) != 1:
        raise ValueError(f"The month spans met versions {list(versions)}: not supported yet")
    met_folder = MET_ROOT / versions[0]
    write_met_restore_job(
        month_folder / "restore_met.sbatch",
        month_folder,
        versions[0],
        prepared["met_start"].min(),
        prepared["met_end"].max(),
        met_folder,
    )
    lines = [
        "#!/bin/bash",
        "# Restore the met; then for each day the surface-pressure run, then the main run.",
        "# afterok: a job starts only if the one before it succeeded.",
        "set -e",
        f"cd {month_folder}",
        "restore=$(sbatch --parsable restore_met.sbatch)",
        'echo "met restore: job $restore"',
    ]
    for day in prepared["day"]:
        folder = pd.Timestamp(day).strftime("%Y%m%d")
        lines += [
            f"ps=$(sbatch --parsable --dependency=afterok:$restore "
            f"{folder}/surface_pressure/job.sbatch)",
            f"main=$(sbatch --parsable --dependency=afterok:$ps {folder}/job.sbatch)",
            f'echo "{folder}: surface pressure job $ps, main job $main"',
        ]
    (month_folder / "submit.sh").write_text("\n".join(lines) + "\n")
    print(summary.to_string(index=False))
    print(f"Prepared {len(prepared)} days in {month_folder}; submit with: bash submit.sh")
    return summary


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if args and args[0] == "main-input":
        write_main_input(args[1])
    else:
        days = [a.removeprefix("--days=") for a in sys.argv[1:] if a.startswith("--days=")]
        prepare_month(*args[:4], test="--test" in sys.argv[1:],
                      days=days[0].split(",") if days else None)  # fmt: skip
