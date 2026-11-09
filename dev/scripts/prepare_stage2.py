"""Prepare stage 2 of the NAME 8.7.1 validation: 1-7 Jan 2013 GOSAT-BRAZIL.

Two purposes, sharing the same met (2 Dec 2012 - 7 Jan 2013):

    day_<YYYYMMDD>  (2-7 Jan): one run per day with the legacy 2019 window, like run B of stage 1,
                    to compare with the production footprints on more days, including the busiest
                    (7 Jan, 64 soundings, 1,088 sources).
    batch_0107:     ONE run for all 7 days (2,720 sources) with the corrected window
                    (``all_releases``), to test batching (batching.py). Compare it with the daily
                    runs (and stage 1 run A for 1 Jan), against the B-B2 noise yardstick.

Releases reproduce the legacy sources exactly, including the widened boxes of binned soundings
(taken from the legacy input files, which is the only place they are recorded), and the legacy
surface pressure (D16). Every plan is checked against the legacy input file before writing.

Also writes a met-restore job; the NAME jobs must wait for it (submit with
--dependency=afterok:<restore job>).

    uv run python dev/scripts/prepare_stage2.py /work/scratch-pw5/mattrigby/name-satellite/stage2_20261109
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from name_satellite.batching import combine_release_plans
from name_satellite.legacy import legacy_surface_pressure_hpa, read_name_input_sources
from name_satellite.lotus import write_job_script
from name_satellite.met import met_version_for_period
from name_satellite.name_input import (
    DispersionDomain,
    OutputGrid,
    RunSettings,
    run_timing,
    write_name_input,
)
from name_satellite.refdata import ref_path
from name_satellite.releases import build_release_plan

SCRATCH = Path("/work/scratch-pw5/mattrigby/name-satellite")
REPO = Path(__file__).resolve().parents[2]
LEGACY_RUN = ref_path("conc2018_input").parents[1]  # all 31 days; sources identical to 2019 run
RUN_NAME = "SOUTHAMERICA_GOSAT-BRAZIL_column"
DAYS = [f"201301{d:02d}" for d in range(1, 8)]
LEVELS = 17

RESTORE_JOB = """\
#!/bin/bash
#SBATCH --job-name=restore_met
#SBATCH --account=acrg
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=08:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output={work}/restore-%j.out
#SBATCH --error={work}/restore-%j.err

cd {repo}
time ~/.local/bin/uv run --no-sync python dev/scripts/restore_met.py {version} {start} {end} {folder}
"""


def legacy_day_plan(day):
    """Release plan for one day that reproduces the legacy sources; checked against them."""
    obs = xr.open_dataset(LEGACY_RUN / "Observations" / f"gosat-fts_gosat_{day}_ch4-column.nc")
    legacy = read_name_input_sources(LEGACY_RUN / "Input_files" / f"BackRun_{RUN_NAME}_{day}.txt")
    boxes = (legacy["dx"].to_numpy()[::LEVELS], legacy["dy"].to_numpy()[::LEVELS])
    plan = build_release_plan(
        obs, surface_pressure_hpa=legacy_surface_pressure_hpa(legacy), box_degrees=boxes
    )
    # Fail loudly if anything differs from the legacy sources. The legacy input file has ~7
    # significant figures: pressures agree to within 0.003 Pa, positions to ~1e-7 relative.
    if plan["source_name"].tolist() != legacy["name"].tolist():
        raise ValueError(f"{day}: source names differ from legacy")
    if not (plan["release_end"] == legacy["stop_time"]).all():
        raise ValueError(f"{day}: release times differ from legacy")
    for ours, theirs in [("lon", "x"), ("lat", "y"), ("dlon", "dx"), ("dlat", "dy")]:
        if not np.allclose(plan[ours], legacy[theirs], rtol=1e-6, atol=0):
            raise ValueError(f"{day}: {ours} differs from legacy {theirs}")
    for ours, theirs in [("z_pa", "z"), ("dz_pa", "dz")]:
        if not np.allclose(plan[ours], legacy[theirs], rtol=0, atol=0.01):  # Pa
            raise ValueError(f"{day}: {ours} differs from legacy {theirs} by > 0.01 Pa")
    if plan["particles_per_second"].tolist() != legacy["particles_per_second"].tolist():
        raise ValueError(f"{day}: particle rates differ from legacy")
    return plan


def settings_for(label, run, met_version, met_folder, integration, threads):
    return RunSettings(
        run_name=f"stage2_{label}",
        output_folder=str(run / "output"),
        grid=OutputGrid(nx=190, ny=357, dx=0.352, dy=0.234, x_min=-91.329, y_min=-60.981),
        domain=DispersionDomain(x_min=-96.329, x_max=-19.449, y_min=-64.981, y_max=27.557),
        met_version=met_version,
        met_folder=str(met_folder),
        threads=threads,
        integration=integration,
    )


def main(work):
    work = Path(work)
    if work.exists():
        raise SystemExit(f"{work} already exists; choose a new folder")

    plans = {day: legacy_day_plan(day) for day in DAYS}
    met_version = met_version_for_period("2012-12-02", "2013-01-07")
    met_folder = SCRATCH / "met" / met_version.name
    work.mkdir(parents=True)

    runs = []
    # Daily runs, legacy window (1 Jan is stage 1 run B).
    for day in DAYS[1:]:
        run = work / f"day_{day}"
        (run / "output").mkdir(parents=True)
        settings = settings_for(day, run, met_version, met_folder, "legacy_2019", threads=4)
        write_name_input(run / "input.txt", plans[day], settings)
        write_job_script(run / "job.sbatch", run, run / "input.txt", threads=4, qos="high",
                         time="12:00:00", memory="32G", job_name=f"s2_{day}")  # fmt: skip
        runs.append(run)
        print(f"{run.name}: {len(plans[day])} sources")

    # One batched run for all days, corrected window.
    combined = combine_release_plans(plans)
    run = work / "batch_0107"
    (run / "output").mkdir(parents=True)
    settings = settings_for("batch_0107", run, met_version, met_folder, "all_releases", threads=8)
    write_name_input(run / "input.txt", combined, settings)
    combined.to_csv(run / "release_plan.csv", index=False)  # needed to split the output by day
    write_job_script(run / "job.sbatch", run, run / "input.txt", threads=8, qos="high",
                     time="36:00:00", memory="64G", job_name="s2_batch")  # fmt: skip
    runs.append(run)
    timing = run_timing(combined, settings)
    print(f"{run.name}: {len(combined)} sources, window {timing['window_hours']} h")

    # Met for everything any run can reach.
    start = (combined["release_end"].min() - pd.Timedelta(hours=721)).strftime("%Y-%m-%dT%H:%M")
    end = combined["release_start"].max().strftime("%Y-%m-%dT%H:%M")
    (work / "restore_met.sbatch").write_text(
        RESTORE_JOB.format(
            work=work, repo=REPO, version=met_version.name, start=start, end=end, folder=met_folder
        )
    )
    (work / "runs.txt").write_text("\n".join(str(r) for r in runs) + "\n")
    print(f"Prepared {work}; met {start} to {end}")


if __name__ == "__main__":
    main(sys.argv[1])
