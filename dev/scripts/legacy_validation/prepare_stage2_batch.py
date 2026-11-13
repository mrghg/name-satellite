"""Prepare the stage 2 BATCHED run again (1-7 Jan 2013 in one NAME run, 2,720 releases).

The first attempt (stage2_20261109/batch_0107, job 66323622) ran for 8 h and then failed when
writing its fields: all 2,720 releases were in one output group, and NAME allows at most 60,000
characters per output line (about 2,218 releases). The input writer now splits the field requests
into output groups of at most 2,000 (name_input.fields_output_group).

    uv run python dev/scripts/legacy_validation/prepare_stage2_batch.py RUN_FOLDER [--test]

--test: particles travel 1 h only, debug queue (minutes): checks that the split output is written
and reads back, before spending 8 h on the full run. The met is the shared, already restored
Mk6 folder used by prepare_stage2.py.
"""

import dataclasses
import sys
from pathlib import Path

from prepare_stage2 import DAYS, SCRATCH, legacy_day_plan, settings_for

from name_satellite.batching import combine_release_plans
from name_satellite.lotus import write_job_script
from name_satellite.met import met_version_for_period
from name_satellite.name_input import run_timing, write_name_input


def main(run, test=False):
    run = Path(run)
    if run.exists():
        raise SystemExit(f"{run} already exists; choose a new folder")
    combined = combine_release_plans({day: legacy_day_plan(day) for day in DAYS})
    met_version = met_version_for_period("2012-12-02", "2013-01-07")
    met_folder = SCRATCH / "met" / met_version.name

    (run / "output").mkdir(parents=True)
    threads = 8
    settings = settings_for(run.name, run, met_version, met_folder, "all_releases", threads)
    if test:
        settings = dataclasses.replace(settings, max_age_hours=1)
    write_name_input(run / "input.txt", combined, settings)
    combined.to_csv(run / "release_plan.csv", index=False)  # needed to split the output by day
    if test:
        write_job_script(run / "job.sbatch", run, run / "input.txt", threads=threads,
                         qos="debug", time="01:00:00", memory="32G", job_name="s2_btest")  # fmt: skip
    else:
        write_job_script(run / "job.sbatch", run, run / "input.txt", threads=threads,
                         qos="high", time="36:00:00", memory="64G", job_name="s2_batch")  # fmt: skip
    timing = run_timing(combined, settings)
    print(f"{run}: {len(combined)} sources, window {timing['window_hours']} h, test={test}")


if __name__ == "__main__":
    main(sys.argv[1], test="--test" in sys.argv[2:])
