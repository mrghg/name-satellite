"""Write the run record (spec §8) for the main block runs and every other finished run.

    uv run python dev/scripts/glide_tests/write_main_record.py MAIN_FOLDER [OTHER_FOLDER ...]

Reads the run folders and job ids from ``jobs.txt`` in each folder (lines "<run folder> <job>";
runs that failed at start-up and were rerun appear under their ``_try2`` folder). Writes
``MAIN_FOLDER/run_record_NAME_main.yaml`` (restricted: contains NAME's input file) and
``MAIN_FOLDER/lamont_surface_pressure.csv`` (both Lamont blocks).
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd

from name_satellite.glide_runs.record import build_record, write_record

DEPARTURES = [
    "Particles per release: NAME releases 21,601 (6/s x 3,600 + 1) per point release and "
    "100,801 (28/s x 3,600 + 1) per Lamont release, within 'a few particles' of the spec (§4).",
    "Surface type at the receptor (§3, §8): not available; NAME reports Land Fraction as missing "
    "with Mk10 met. The model ground height is given.",
    "Lamont column (§5): one NAME source per release, spread evenly in pressure (measured, job "
    "2296900), not 28 stacked layers.",
    "Lamont blocks: each block was run as two runs of 168 releases (days 1-7 and 8-14) to keep "
    "each job well inside the 48 h limit; settings identical; the delivered file joins them.",
    "Single-release run file names: the spec gives no run name; 'case<n>' (n = §9 standard case) "
    "is used, with '_nomeso' for the nomeso runs.",
    "Six seed runs (MHD seed1-5, CBW-20 seed2) crashed at start-up on one LOTUS2 node "
    "(host1036) before reading their input; they were rerun unchanged on other nodes (_try2).",
]


def finished_runs(folder):
    """(run folder, job id, purpose) for every completed NAME run listed in folder/jobs.txt."""
    runs = []
    for line in (Path(folder) / "jobs.txt").read_text().splitlines():
        if line.startswith("#") or len(line.split()) != 2:
            continue
        name, job = line.split()
        run = Path(folder) / name
        if not (run / "run.json").exists():
            continue  # met restores, delivery and check jobs, surface-pressure folders
        state = subprocess.run(
            ["sacct", "-n", "-X", "-j", job, "-o", "State"], capture_output=True, text=True
        ).stdout.split()
        if state[:1] == ["COMPLETED"]:
            runs.append({"folder": run, "job": job, "purpose": name})
    return runs


if __name__ == "__main__":
    main_folder = Path(sys.argv[1])
    runs = [run for folder in sys.argv[1:] for run in finished_runs(folder)]
    pressures = pd.concat(
        pd.read_csv(main_folder / f"SGP_NAM-{m}_surface_pressure" / "surface_pressure.csv")
        for m in ("Jan", "Jul")
    )
    csv = main_folder / "lamont_surface_pressure.csv"
    pressures.to_csv(csv, index=False)
    notes = {
        "check10": "One-day check at Mace Head: folder check10-2026-11-16 (all checks pass "
        "under spec e934448; check 8 tolerance one unit in the last printed digit).",
        "check6": "Two identical sources in one Fixed run differ (job 2296902), so main uses "
        "Fixed (§9).",
        "case_times": "Standard cases 3-5 (§9) use the release times of GLIDE PR 22 (17 Nov "
        "2026): CBW-20 7 Jan 2019 02:00, LEF 13 Jul 2019 20:00, TAC 13 Jan 2019 15:00.",
    }
    record = build_record(runs, DEPARTURES, notes, surface_pressure_csv=csv)
    path = write_record(main_folder / "run_record_NAME_main.yaml", record)
    print(f"Wrote {path}: {len(runs)} runs")
