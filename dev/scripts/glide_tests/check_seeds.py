"""GLIDE spec §9: confirm that the five repeat-seed runs of a receptor differ from one another.

NAME's `Random` seed cannot be set or reported, so the spec asks for this check. For each
receptor, compares every pair of delivered seed files (seed1 ... seed5): they must not be
identical, and prints the median over releases of the relative difference of footprint totals.

    uv run python dev/scripts/glide_tests/check_seeds.py DELIVERED_FOLDER MHD JFJ CBW-20
"""

import itertools
import sys
from pathlib import Path

import numpy as np
import xarray as xr


def check(folder, receptor):
    paths = sorted(Path(folder).glob(f"{receptor}_NAME_*_seed[1-5].nc"))
    if len(paths) != 5:
        print(f"SKIP {receptor}: {len(paths)} seed files")
        return
    srr = [xr.open_dataset(p)["srr"].values for p in paths]
    totals = [s.sum(axis=(1, 2)) for s in srr]
    identical = [(i + 1, j + 1) for i, j in itertools.combinations(range(5), 2)
                 if np.array_equal(srr[i], srr[j])]  # fmt: skip
    spread = np.median(np.std(totals, axis=0) / np.mean(totals, axis=0))
    print(f"{'FAIL' if identical else 'PASS'} {receptor}: identical pairs {identical or 'none'}; "
          f"median spread of release totals across seeds {spread:.1%}")  # fmt: skip


if __name__ == "__main__":
    for code in sys.argv[2:]:
        check(sys.argv[1], code)
