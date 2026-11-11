"""Small NAME 8.7.1 run to study the particle-fate output (X20, decision D13/D14).

Two releases of the first 1 Jan 2013 sounding (levels 1 and 9, 1,801 particles each), 30-day
travel, fate output in m agl WITH the met columns (``RunSettings.fate_met``), which include the
pressure at each particle's final position.

Questions:
  * X20: does this small run also write more fate records than particles released (as the 30-day
    stage 1 runs did, +2.2%)? If so, it is the minimal example for the NAME team.
  * D14: is the pressure column present and sensible for every record (lateral exits, top exits,
    particles still inside at the end)?

    uv run python dev/scripts/prepare_fate_test.py WORK_FOLDER [--patched]

``--patched``: run the ACRG-patched NAME (``met.PATCHED_NAME_DIR``), whose fate records also
carry each particle's per-source index and fate code, so duplicates can be identified exactly.
"""

import sys
from pathlib import Path

import xarray as xr

from name_satellite.legacy import legacy_surface_pressure_hpa, read_name_input_sources
from name_satellite.lotus import write_job_script
from name_satellite.met import NAME_DIR, PATCHED_NAME_DIR, met_version_for_period
from name_satellite.name_input import DispersionDomain, OutputGrid, RunSettings, write_name_input
from name_satellite.refdata import ref_path
from name_satellite.releases import build_release_plan

SCRATCH = Path("/work/scratch-pw5/mattrigby/name-satellite")
LEGACY_RUN = ref_path("conc2018_input").parents[1]
RUN_NAME = "SOUTHAMERICA_GOSAT-BRAZIL_column"
DAY = "20130101"
LEVELS = (1, 9)


def main(work, patched=False):
    work = Path(work)
    if work.exists():
        raise SystemExit(f"{work} already exists; choose a new folder")
    obs = xr.open_dataset(LEGACY_RUN / "Observations" / f"gosat-fts_gosat_{DAY}_ch4-column.nc")
    legacy = read_name_input_sources(LEGACY_RUN / "Input_files" / f"BackRun_{RUN_NAME}_{DAY}.txt")
    plan = build_release_plan(obs, surface_pressure_hpa=legacy_surface_pressure_hpa(legacy))
    plan = plan[(plan["point_number"] == 1) & plan["level_number"].isin(LEVELS)].copy()
    plan["source_name"] = [f"SourceID_{n}" for n in range(1, len(plan) + 1)]

    met_version = met_version_for_period("2012-12-02", "2013-01-01")
    settings = RunSettings(
        run_name="fate_test",
        output_folder=str(work / "output"),
        grid=OutputGrid(nx=190, ny=357, dx=0.352, dy=0.234, x_min=-91.329, y_min=-60.981),
        domain=DispersionDomain(x_min=-96.329, x_max=-19.449, y_min=-64.981, y_max=27.557),
        met_version=met_version,
        met_folder=str(SCRATCH / "met" / met_version.name),
        threads=1,
        release_met=False,
        fate_met=True,
    )
    (work / "output").mkdir(parents=True)
    write_name_input(work / "input.txt", plan, settings)
    plan.to_csv(work / "release_plan.csv", index=False)
    write_job_script(work / "job.sbatch", work, work / "input.txt", threads=1, qos="high",
                     time="04:00:00", memory="16G", job_name="fate_test",
                     name_dir=PATCHED_NAME_DIR if patched else NAME_DIR)  # fmt: skip
    # NAME releases rate x 60 s + 1 particles per source (legacy files: exactly 60 x rate + 1).
    expected = int((plan["particles_per_second"] * 60 + 1).sum())
    print(f"Prepared {work}: {len(plan)} releases, {expected} particles expected")


if __name__ == "__main__":
    main(sys.argv[1], patched="--patched" in sys.argv[2:])
