"""Test of release-point met from NAME 8.7.1 (release_met.py) against the legacy met run.

A short run (debug queue): two dummy releases (level 17 of the first and last soundings, 1 h
travel) so that the run covers every sounding time, plus the release-met requests for every
sounding of one January 2013 day. To compare with the legacy met files at the two levels that
can be underground (X9), there are two groups: group 1 samples the level-1 release pressures, group 2 (a copy of the
same soundings) the level-2 pressures.

    uv run python dev/scripts/prepare_release_met_test.py WORK_FOLDER [YYYYMMDD]   (default 20130101)

1 Jan has no underground releases; 2 Jan has three (level 1 of soundings 23, 25, 26).
"""

import sys
from pathlib import Path

import pandas as pd
import xarray as xr

from name_satellite.legacy import legacy_surface_pressure_hpa, read_name_input_sources
from name_satellite.lotus import write_job_script
from name_satellite.met import met_version_for_period
from name_satellite.name_input import DispersionDomain, OutputGrid, RunSettings, name_input_text
from name_satellite.refdata import ref_path
from name_satellite.release_met import release_met_blocks, release_met_soundings
from name_satellite.releases import build_release_plan

SCRATCH = Path("/work/scratch-pw5/mattrigby/name-satellite")
LEGACY_RUN = ref_path("conc2018_input").parents[1]  # all January days
RUN_NAME = "SOUTHAMERICA_GOSAT-BRAZIL_column"


def main(work, day="20130101"):
    work = Path(work)
    if work.exists():
        raise SystemExit(f"{work} already exists; choose a new folder")

    obs = xr.open_dataset(LEGACY_RUN / "Observations" / f"gosat-fts_gosat_{day}_ch4-column.nc")
    legacy = read_name_input_sources(LEGACY_RUN / "Input_files" / f"BackRun_{RUN_NAME}_{day}.txt")
    plan = build_release_plan(obs, surface_pressure_hpa=legacy_surface_pressure_hpa(legacy))

    level1 = release_met_soundings(plan)
    level2_plan = plan[plan["level_number"] == 2].assign(level_number=1)
    level2 = release_met_soundings(level2_plan)
    level2["sounding"] = level2["sounding"].str.replace("S", "L2S")
    level2["group"] = 2
    soundings = pd.concat([level1, level2], ignore_index=True)

    # Dummy releases: top level of the first and last soundings, travelling 1 h.
    dummy = plan[
        (plan["level_number"] == 17) & plan["point_number"].isin([1, plan["point_number"].max()])
    ]
    dummy = dummy.assign(source_name=["SourceID_1", "SourceID_2"])

    met_version = met_version_for_period(day, day)
    settings = RunSettings(
        run_name=f"release_met_test_{day}",
        output_folder=str(work / "output"),
        grid=OutputGrid(nx=190, ny=357, dx=0.352, dy=0.234, x_min=-91.329, y_min=-60.981),
        domain=DispersionDomain(x_min=-96.329, x_max=-19.449, y_min=-64.981, y_max=27.557),
        met_version=met_version,
        met_folder=str(SCRATCH / "met" / met_version.name),
        max_age_hours=1,
    )
    (work / "output").mkdir(parents=True)
    text = name_input_text(dummy, settings) + "\n" + release_met_blocks(soundings, 0.352, 0.234)
    (work / "input.txt").write_text(text)
    soundings.to_csv(work / "soundings.csv", index=False)
    write_job_script(work / "job.sbatch", work, work / "input.txt", threads=1, qos="debug",
                     time="00:30:00", memory="8G", job_name="release_met_test")  # fmt: skip
    print(f"Prepared {work}: {len(soundings)} sampled soundings")


if __name__ == "__main__":
    main(*sys.argv[1:])
