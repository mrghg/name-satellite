"""Prepare stage 1 of the NAME 8.7.1 validation: 1 Jan 2013 GOSAT-BRAZIL, two runs.

    run_A: integration window covering every release (our default, X18 fixed)
    run_B: the legacy 2019 integration window (721 h from the first sounding)

Everything else is identical: all 425 releases, legacy surface pressure (D16), 30-day back-run,
legacy domain/grid, Mk6 met, fixed random seed, 4 threads on the `high` QoS. Because the seed is
fixed and the thread count doesn't change results, A - B is purely the effect of the window.

Also writes a met-restore job (run first; the NAME jobs depend on it):

    <work>/restore_met.sbatch
    <work>/run_A/{input.txt, job.sbatch, output/}
    <work>/run_B/{...}

    uv run python dev/scripts/prepare_stage1.py /work/scratch-pw5/mattrigby/name-satellite/stage1_20261109
"""

import sys
from pathlib import Path

import pandas as pd
import xarray as xr

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


def main(work):
    work = Path(work)
    if work.exists():
        raise SystemExit(f"{work} already exists; choose a new folder")
    work.mkdir(parents=True)

    obs = xr.open_dataset(ref_path("obs"))
    legacy_sources = read_name_input_sources(ref_path("mr2019_input"))
    plan = build_release_plan(obs, surface_pressure_hpa=legacy_surface_pressure_hpa(legacy_sources))
    met_version = met_version_for_period("2012-12-02", "2013-01-01")
    met_folder = SCRATCH / "met" / met_version.name

    for label, integration in (("A", "all_releases"), ("B", "legacy_2019")):
        run = work / f"run_{label}"
        (run / "output").mkdir(parents=True)
        settings = RunSettings(
            run_name=f"stage1_20130101_{label}",
            output_folder=str(run / "output"),
            grid=OutputGrid(nx=190, ny=357, dx=0.352, dy=0.234, x_min=-91.329, y_min=-60.981),
            domain=DispersionDomain(x_min=-96.329, x_max=-19.449, y_min=-64.981, y_max=27.557),
            met_version=met_version,
            met_folder=str(met_folder),
            threads=4,
            integration=integration,
        )
        write_name_input(run / "input.txt", plan, settings)
        write_job_script(
            run / "job.sbatch",
            run,
            run / "input.txt",
            threads=4,
            qos="high",
            time="12:00:00",
            memory="16G",
            job_name=f"stage1_{label}",
        )
        timing = run_timing(plan, settings)
        print(f"run_{label}: window from {timing['window_start']}, {timing['window_hours']} h")

    # Met for the whole period both runs can reach.
    start = (plan["release_end"].min() - pd.Timedelta(hours=720)).strftime("%Y-%m-%dT%H:%M")
    end = plan["release_start"].max().strftime("%Y-%m-%dT%H:%M")
    (work / "restore_met.sbatch").write_text(
        RESTORE_JOB.format(
            work=work, repo=REPO, version=met_version.name, start=start, end=end, folder=met_folder
        )
    )
    print(f"Prepared {work}; met {start} to {end}")


if __name__ == "__main__":
    main(sys.argv[1])
