"""Stage 1 noise test: run B again with a different random seed (run_B2).

Run B (prepare_stage1.py) used ``Random Seed = Fixed``. B2 is identical except for
``Random (Parallel)``, the legacy setting, so B - B2 is pure particle noise. If NAME 8.7.1
reproduces legacy, B - production should be no larger than B - B2.

Writes <work>/run_B2/{input.txt, job.sbatch, output/} into an existing stage 1 folder, after
checking that every met file the run needs is already restored (no restore job this time).

    uv run python dev/scripts/legacy_validation/prepare_stage1_noise.py /work/scratch-pw5/mattrigby/name-satellite/stage1_20261109
"""

import sys
from pathlib import Path

import pandas as pd
import xarray as xr

from name_satellite.lotus import write_job_script
from name_satellite.met import met_version_for_period, required_met_files
from name_satellite.name_input import (
    DispersionDomain,
    OutputGrid,
    RunSettings,
    run_timing,
    write_name_input,
)
from name_satellite.releases import build_release_plan
from name_satellite.validation.legacy_files import (
    legacy_surface_pressure_hpa,
    read_name_input_sources,
)
from name_satellite.validation.refdata import ref_path

SCRATCH = Path("/work/scratch-pw5/mattrigby/name-satellite")


def main(work):
    work = Path(work)
    run = work / "run_B2"
    if run.exists():
        raise SystemExit(f"{run} already exists; choose a new folder")

    obs = xr.open_dataset(ref_path("obs"))
    legacy_sources = read_name_input_sources(ref_path("mr2019_input"))
    plan = build_release_plan(obs, surface_pressure_hpa=legacy_surface_pressure_hpa(legacy_sources))
    met_version = met_version_for_period("2012-12-02", "2013-01-01")
    met_folder = SCRATCH / "met" / met_version.name

    # Same settings as run_B in prepare_stage1.py, except the random seed.
    settings = RunSettings(
        run_name="stage1_20130101_B2",
        output_folder=str(run / "output"),
        grid=OutputGrid(nx=190, ny=357, dx=0.352, dy=0.234, x_min=-91.329, y_min=-60.981),
        domain=DispersionDomain(x_min=-96.329, x_max=-19.449, y_min=-64.981, y_max=27.557),
        met_version=met_version,
        met_folder=str(met_folder),
        threads=4,
        integration="legacy_2019",
        random_seed="Random (Parallel)",
    )

    # The met must already be in scratch (restored for runs A and B).
    timing = run_timing(plan, settings)
    window_end = timing["window_start"] + pd.Timedelta(hours=timing["window_hours"])
    needed = required_met_files(met_version, timing["window_start"], window_end)
    missing = [name for name in needed if not (met_folder / name).exists()]
    if missing:
        raise SystemExit(f"{len(missing)} met files not restored, e.g. {missing[0]}")

    (run / "output").mkdir(parents=True)
    write_name_input(run / "input.txt", plan, settings)
    write_job_script(
        run / "job.sbatch",
        run,
        run / "input.txt",
        threads=4,
        qos="high",
        time="12:00:00",
        memory="16G",
        job_name="stage1_B2",
    )
    print(f"Prepared {run}: window from {timing['window_start']}, {timing['window_hours']} h")


if __name__ == "__main__":
    main(sys.argv[1])
