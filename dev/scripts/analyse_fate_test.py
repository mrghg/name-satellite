"""Analyse the particle-fate output of the small fate test run (prepare_fate_test.py) for X20.

Prints: records vs particles released, the fate of each record (lateral exit / top exit / end of
run), the 10-second near-duplicate pairs, and when (clock time) the records were written.

    uv run python dev/scripts/analyse_fate_test.py RUN_FOLDER
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from name_satellite.name_output import read_name_fate_file

# Dispersion domain of the run (prepare_fate_test.py) and its top (m agl).
X_MIN, X_MAX, Y_MIN, Y_MAX, Z_TOP = -96.329, -19.449, -64.981, 27.557, 19000.0
EDGE_TOLERANCE = 0.002  # degrees: exits are recorded up to ~0.001 deg inside an edge (rounding)
MAX_AGE_HOURS = 720.0
RELEASE_MIDDLE_HOURS = 14 + 59 / 60 + 20 / 3600  # releases from 14:58:50 to 14:59:50 (UTC)


def main(run):
    run = Path(run)
    plan = pd.read_csv(run / "release_plan.csv")
    released = int((plan["particles_per_second"] * 60 + 1).sum())
    fates, z_units = read_name_fate_file(run / "output" / "Particles_C1.txt")
    assert z_units == "m agl"
    print(f"records {len(fates)}, particles released {released}, "
          f"excess {len(fates) - released} ({100 * (len(fates) / released - 1):.1f}%)")  # fmt: skip
    print("records per source:", fates.groupby("source_number").size().to_dict(),
          "expected", (plan["particles_per_second"] * 60 + 1).tolist())  # fmt: skip

    lon, lat, z, age = fates["lon"], fates["lat"], fates["z"], fates["age_hours"]
    lateral = (
        (lon <= X_MIN + EDGE_TOLERANCE) | (lon >= X_MAX - EDGE_TOLERANCE)
        | (lat <= Y_MIN + EDGE_TOLERANCE) | (lat >= Y_MAX - EDGE_TOLERANCE)
    )  # fmt: skip
    top = ~lateral & (z >= Z_TOP - 1)
    end = ~lateral & ~top & (age >= MAX_AGE_HOURS - 1)
    other = ~(lateral | top | end)
    print(f"fates: lateral {lateral.sum()}, top {top.sum()}, end of run {end.sum()}, "
          f"other {other.sum()}")  # fmt: skip

    # Near-duplicates: same source, one 10-s time step apart, almost the same position.
    pairs = 0
    for _, group in fates.groupby("source_number"):
        values = group[["age_hours", "lon", "lat", "z"]].to_numpy()
        for i in range(len(values)):
            d = np.abs(values - values[i])
            close = (d[:, 0] > 0) & (d[:, 0] < 0.2) & (d[:, 1] < 0.01) & (d[:, 2] < 0.01)
            close &= d[:, 3] < 5
            pairs += int(close[i + 1 :].sum())
    print(f"near-duplicate pairs (same source, < 0.2 h, < 0.01 deg, < 5 m): {pairs}")

    # Clock time at which each record was written (release time - travel time, +-30 s).
    clock_minute = ((RELEASE_MIDDLE_HOURS - age) * 60) % 60
    counts = np.histogram(clock_minute, bins=30, range=(0, 60))[0]
    # Records are written only at some phases of the 6-minute sync cycle (empty bins every 6 min);
    # the excess shows as larger counts in the bins just before each whole hour.
    print("records by clock minute of the hour (2-minute bins, 00-02 ... 58-60):")
    print(" ", counts.tolist())


if __name__ == "__main__":
    main(sys.argv[1])
