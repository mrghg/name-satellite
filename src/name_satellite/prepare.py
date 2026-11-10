"""Prepare a month of NAME runs from Level-2 satellite data: workflow steps 1-4 chained together.

Where this sits (see docs/workflow.md)
--------------------------------------
    L2 files (CEDA)  ->  **prepare**: per day, observation file + release plan + NAME input + job
    ->  run the jobs on LOTUS2  ->  read and process the output (steps 5-6)

For each day of the month, ``prepare_day``:
  1. downloads the day's L2 file (``l2_files``; days without a file are skipped and reported);
  2. reads it and keeps good soundings in the site's box (``observations``: quality flag, defined
     pressure levels, pressures decreasing from the surface, lat/lon bounds), binned if the site
     says so;
  3. writes the observation file, keyed by sounding ID (D18);
  4. plans the releases (``releases.build_release_plan``): bottom layer from the RETRIEVAL's
     surface pressure (the legacy separate surface-pressure run does not exist for new data; the
     run's release-point met reports the model's surface pressure, so X9 can be checked);
  5. writes the NAME input (with release-point met, X21) and the Slurm job.

``prepare_month`` also writes one met-restore job for the whole month and ``submit.sh``, which
submits the restore job and then every day's NAME job after it (Slurm ``afterok``).

Run folder layout (``WORK/<site>_<product>_<YYYYMM>/``)::

    restore_met.sbatch, submit.sh, summary.csv
    YYYYMMDD/obs_<site>_<YYYYMMDD>.nc, release_plan.csv, input.txt, job.sbatch, output/

    uv run python -m name_satellite.prepare "GOSAT-1 v10.0" GOSAT-BRAZIL 2013-01 WORK [--test]

``--test``: particles travel 1 hour only, debug queue, to check a set-up quickly (lesson of X20).
"""

import dataclasses
import sys
from pathlib import Path

import pandas as pd

from name_satellite.domains import SITES
from name_satellite.l2_files import L2FileMissing, download_l2
from name_satellite.lotus import write_job_script, write_met_restore_job
from name_satellite.met import met_version_for_period
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
from name_satellite.release_met import MAX_SOUNDINGS
from name_satellite.releases import build_release_plan

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
    )
    if test:
        settings = dataclasses.replace(settings, max_age_hours=1)
    return settings


def prepare_day(product, site, day, month_folder, threads=THREADS, test=False):
    """Prepare one day's NAME run. Returns a summary dict (one row of summary.csv)."""
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
    obs = to_obs(table)
    boxes = (obs["dlon"].values, obs["dlat"].values) if "dlat" in obs else None
    plan = build_release_plan(obs, footprint_diameter_km=site.footprint_km, box_degrees=boxes)
    plan.to_csv(run_folder / "release_plan.csv", index=False)

    settings = run_settings(site, run_folder, day, threads, test)
    write_name_input(run_folder / "input.txt", plan, settings)
    qos, time = ("debug", "00:30:00") if test else ("high", "12:00:00")
    write_job_script(run_folder / "job.sbatch", run_folder, run_folder / "input.txt",
                     threads=threads, qos=qos, time=time, memory="32G",
                     job_name=f"{site.name}_{day:%Y%m%d}")  # fmt: skip
    timing = run_timing(plan, settings)
    row.update(
        status="prepared",
        sources=len(plan),
        met_version=settings.met_version.name,
        met_start=f"{timing['window_start']:%Y-%m-%dT%H:%M}",
        met_end=f"{plan['release_start'].max():%Y-%m-%dT%H:%M}",
    )
    return row


def prepare_month(product, site_name, month, work, threads=THREADS, test=False):
    """Prepare every day of ``month`` ("YYYY-MM"). Returns the summary table."""
    site = SITES[site_name]
    month = pd.Period(month, freq="M")
    month_folder = Path(work) / f"{site.name}_{product.replace(' ', '_')}_{month:%Y%m}"
    if month_folder.exists():
        raise SystemExit(f"{month_folder} already exists; choose a new folder")
    month_folder.mkdir(parents=True)

    rows = [
        prepare_day(product, site, day, month_folder, threads, test)
        for day in pd.date_range(month.start_time, month.end_time.normalize(), freq="D")
    ]
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
        "# Restore the met, then each day's NAME job after it (afterok: if the restore succeeded).",
        "set -e",
        f"cd {month_folder}",
        "restore=$(sbatch --parsable restore_met.sbatch)",
        'echo "met restore: job $restore"',
    ]
    for day in prepared["day"]:
        folder = pd.Timestamp(day).strftime("%Y%m%d")
        lines.append(f'echo "{folder}: job $(sbatch --parsable --dependency=afterok:$restore '
                     f'{folder}/job.sbatch)"')  # fmt: skip
    (month_folder / "submit.sh").write_text("\n".join(lines) + "\n")
    print(summary.to_string(index=False))
    print(f"Prepared {len(prepared)} days in {month_folder}; submit with: bash submit.sh")
    return summary


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    prepare_month(*args[:4], test="--test" in sys.argv[1:])
