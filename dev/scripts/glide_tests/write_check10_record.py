"""Write the run record (spec §8) for the §10 one-day check at Mace Head, 16 Nov 2026.

    uv run python dev/scripts/glide_tests/write_check10_record.py CHECK_FOLDER

Reads ``CHECK_FOLDER/check10_results.txt`` (from ``check10.py``) and the four run folders; writes
``CHECK_FOLDER/run_record_NAME_check10.yaml`` (restricted: contains NAME's input file).
"""

import sys
from pathlib import Path

from name_satellite.glide_runs.record import build_record, write_record

RUNS = [
    ("A-fixed", 2306225, "§10 one-day check; delivered as main"),
    ("B-fixed-repeat", 2306226, "§10 check 5: repeat of A (Fixed)"),
    ("C-random", 2306227, "§10 check 5: Random seed"),
    ("D-fixed-unpatched", 2306228, "§10 check 8: stock NAME 8.7.1"),
]

DEPARTURES = [
    "Particles per release: NAME releases 21,601 (6/s x 3,600 + 1), within 'a few particles' "
    "of 21,600 (§4).",
    "Surface type at the receptor (§3, §8): not available; NAME reports Land Fraction as missing "
    "with Mk10 met. The model ground height is given.",
    "§10 check 8: the patched and stock builds give identical particle-fate records, but one "
    "footprint value (of 24 x 44,361) differs in the 7th significant figure (10.50166 vs "
    "10.50167 ppm s). Not bit-identical; believed to be floating-point rounding between two "
    "separately built executables.",
    "Lamont column (§5): one NAME source per release, spread evenly in pressure (measured, job "
    "2296900), not 28 stacked layers.",
    "Single-release run file names: the spec gives no run name; 'case<n>' (n = §9 standard case) "
    "is used until it does.",
]

NOTES = {
    "check6": "Two identical sources in one Fixed run differ (test of 16 Nov 2026, job "
    "2296902): 0 of 21,601 end positions equal. So main uses Fixed (§9).",
    "footprint_size": "Totals (median 3.9 m2 s/mol) are below the 5-17 of the group's earlier "
    "Mace Head files, but not by a factor of ten. Likely causes: the small EUR box (87% of "
    "particles leave through the west edge, 30 W) and the 100 m minimum boundary-layer depth "
    "(earlier files: 40 m). The 'group' run will measure the settings' part.",
    "memory_and_time": "One-day runs (24 releases, 240 h): 21-31 min on 8 threads, about 5 GB. "
    "Memory probe (job 2301599): 5.2 M live particles in 5.1 GB.",
}

if __name__ == "__main__":
    folder = Path(sys.argv[1])
    runs = [{"folder": folder / name, "job": job, "purpose": why} for name, job, why in RUNS]
    notes = {"check10_results": (folder / "check10_results.txt").read_text(), **NOTES}
    path = write_record(
        folder / "run_record_NAME_check10.yaml", build_record(runs, DEPARTURES, notes)
    )
    print(f"Wrote {path}")
