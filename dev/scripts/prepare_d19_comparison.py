"""Measure the effect of D19: release layers scaled to the model surface vs the retrieval's levels.

1 Jan 2013, GOSAT-1 v10.0, GOSAT-BRAZIL, full 30-day runs, same random seed. Two runs from the
same observation file:

  * ``scaled/20130101``: the normal workflow (``prepare.prepare_day``): met-only surface-pressure
    run, then the main job, which scales every level to the model surface (D19);
  * ``retrieval_levels/20130101``: the retrieval's own levels (the set-up before D19), otherwise
    identical.

Releases differ only in their pressures (``z_pa``, ``dz_pa``), so with the same seed the
footprint differences are the effect of D19, plus the noise from particles taking different
paths.

    uv run python dev/scripts/prepare_d19_comparison.py WORK_FOLDER
    bash WORK_FOLDER/submit.sh
    uv run python dev/scripts/compare_d19.py WORK_FOLDER     (when both have finished)

Assumes the met for December 2012 - 1 January 2013 is already restored (it is, from the January
2013 test runs) in ``prepare.MET_ROOT``.
"""

import shutil
import sys
from pathlib import Path

import pandas as pd
import xarray as xr

from name_satellite import prepare
from name_satellite.domains import SITES
from name_satellite.lotus import write_job_script
from name_satellite.name_input import write_name_input
from name_satellite.releases import build_release_plan

PRODUCT = "GOSAT-1 v10.0"
SITE = SITES["GOSAT-BRAZIL"]
DAY = "2013-01-01"


def main(work):
    work = Path(work)
    if work.exists():
        raise SystemExit(f"{work} already exists; choose a new folder")

    # (c) the normal workflow, D19.
    row = prepare.prepare_day(PRODUCT, SITE, DAY, work / "scaled")
    print(row)
    scaled = work / "scaled" / "20130101"

    # (a) the same soundings with the retrieval's own levels.
    folder = work / "retrieval_levels" / "20130101"
    (folder / "output").mkdir(parents=True)
    obs_name = f"obs_{SITE.name}_20130101.nc"
    shutil.copy(scaled / obs_name, folder / obs_name)
    with xr.open_dataset(folder / obs_name) as f:
        obs = f.load()
    plan = build_release_plan(obs, footprint_diameter_km=SITE.footprint_km,
                              box_degrees=prepare.release_boxes(obs))  # fmt: skip
    plan.to_csv(folder / "release_plan.csv", index=False)
    settings = prepare.run_settings(SITE, folder, pd.Timestamp(DAY), prepare.THREADS,
                                    test=False)  # fmt: skip
    write_name_input(folder / "input.txt", plan, settings)
    write_job_script(folder / "job.sbatch", folder, folder / "input.txt", threads=prepare.THREADS,
                     qos="high", time="12:00:00", memory="32G",
                     job_name="d19_retrieval_levels")  # fmt: skip

    lines = [
        "#!/bin/bash",
        "set -e",
        f"cd {work}",
        "ps=$(sbatch --parsable scaled/20130101/surface_pressure/job.sbatch)",
        "main=$(sbatch --parsable --dependency=afterok:$ps scaled/20130101/job.sbatch)",
        "control=$(sbatch --parsable retrieval_levels/20130101/job.sbatch)",
        'echo "surface pressure $ps, scaled (D19) $main, retrieval levels $control"',
    ]
    (work / "submit.sh").write_text("\n".join(lines) + "\n")
    print(f"Prepared in {work}; submit with: bash {work / 'submit.sh'}")


if __name__ == "__main__":
    main(sys.argv[1])
