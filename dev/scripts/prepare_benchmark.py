"""Prepare a NAME 8.7.1 threading benchmark on LOTUS2 (debug queue).

The same run (all 425 releases of 1 Jan 2013, legacy surface pressure (D16), 24 h maximum particle
age, Mk6 met, timers on) with 1, 4 and 8 OpenMP threads, each in its own folder:

    <work>/threads_<n>/{input.txt, job.sbatch, output/}

Compare elapsed time, CPU efficiency (sacct) and NAME's timer summaries to choose a thread count.

    uv run python dev/scripts/prepare_benchmark.py /work/scratch-pw5/mattrigby/name-satellite/benchmark_20261109
    for n in 1 4 8; do sbatch <work>/threads_$n/job.sbatch; done
"""

import sys
from pathlib import Path

import xarray as xr

from name_satellite.legacy import legacy_surface_pressure_hpa, read_name_input_sources
from name_satellite.lotus import write_job_script
from name_satellite.met import met_version_for_period
from name_satellite.name_input import DispersionDomain, OutputGrid, RunSettings, write_name_input
from name_satellite.refdata import ref_path
from name_satellite.releases import build_release_plan

SCRATCH = Path("/work/scratch-pw5/mattrigby/name-satellite")
THREADS = (1, 4, 8)


def main(work):
    work = Path(work)
    if work.exists():
        raise SystemExit(f"{work} already exists; choose a new folder")

    obs = xr.open_dataset(ref_path("obs"))
    legacy_sources = read_name_input_sources(ref_path("mr2019_input"))
    plan = build_release_plan(obs, surface_pressure_hpa=legacy_surface_pressure_hpa(legacy_sources))

    for threads in THREADS:
        run = work / f"threads_{threads}"
        (run / "output").mkdir(parents=True)
        settings = RunSettings(
            run_name=f"benchmark_20130101_threads{threads}",
            output_folder=str(run / "output"),
            grid=OutputGrid(nx=190, ny=357, dx=0.352, dy=0.234, x_min=-91.329, y_min=-60.981),
            domain=DispersionDomain(x_min=-96.329, x_max=-19.449, y_min=-64.981, y_max=27.557),
            met_version=met_version_for_period("2012-12-31", "2013-01-01"),
            met_folder=str(SCRATCH / "met" / "UMG_Mk6PT"),
            max_age_hours=24,
            threads=threads,
        )
        write_name_input(run / "input.txt", plan, settings)
        write_job_script(
            run / "job.sbatch",
            run,
            run / "input.txt",
            threads=threads,
            qos="debug",
            time="01:00:00",
            memory="16G",
            job_name=f"name_bench_{threads}",
        )
    print(f"Prepared {work}")


if __name__ == "__main__":
    main(sys.argv[1])
