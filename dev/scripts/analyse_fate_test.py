"""Analyse the particle-fate output of the small fate test run (prepare_fate_test.py) for X20.

Prints: records vs particles released, the fate of each record (lateral exit / top exit / end of
run), the 10-second near-duplicate pairs, and when (clock time) the records were written.

For a run of the ACRG-patched NAME (``prepare_fate_test.py --patched``), whose records carry each
particle's per-source index and fate code, it also checks that every particle 1..N of each source
appears, lists the particles recorded more than once, compares the fate codes with the fates
inferred from position, and (with ``--stock RUN``) compares the records with a stock run of the
same input: with the same random seed, the particles should be in the same places.

    uv run python dev/scripts/analyse_fate_test.py RUN_FOLDER [--stock STOCK_RUN_FOLDER]
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


def patched_checks(fates, plan, lateral, top, end):
    """Checks that need the patched NAME's particle index and fate code."""
    expected = (plan["particles_per_second"] * 60 + 1).astype(int).tolist()
    for source, n in enumerate(expected, start=1):
        index = fates.loc[fates["source_number"] == source, "particle_index"]
        missing = sorted(set(range(1, n + 1)) - set(index))
        extra = sorted(set(index) - set(range(1, n + 1)))
        repeated = index[index.duplicated()].nunique()
        print(f"source {source}: {n} released, {index.nunique()} distinct particles recorded, "
              f"{repeated} recorded more than once, missing {len(missing)}, outside 1..N "
              f"{len(extra)}")  # fmt: skip

    key = ["source_number", "particle_index"]
    repeats = fates[fates.duplicated(key, keep=False)].sort_values([*key, "age_hours"])
    if len(repeats):
        steps = repeats.groupby(key)["age_hours"].agg(lambda a: (a.max() - a.min()) * 3600)
        print(
            f"repeated particles: {repeats.groupby(key).ngroups}; seconds between first and "
            f"last record: {sorted(steps.round().astype(int).value_counts().to_dict().items())}"
        )
        print("fate codes of repeated records:", repeats["fate_code"].value_counts().to_dict())
        print("first repeats:")
        print(repeats.head(6).to_string(index=False))

    inferred = np.select([lateral, top, end], ["lateral", "top", "end"], "other")
    print("fate code by fate inferred from position:")
    print(pd.crosstab(fates["fate_code"], inferred).to_string())


def compare_with_stock(fates, stock_run):
    """Records of the patched run vs a stock run of the same input (same seed)."""
    stock, _ = read_name_fate_file(Path(stock_run) / "output" / "Particles_C1.txt")
    columns = ["source_number", "age_hours", "lon", "lat", "z"]
    a = fates[columns].sort_values(columns).reset_index(drop=True)
    b = stock[columns].sort_values(columns).reset_index(drop=True)
    same = len(a) == len(b) and bool((a == b).all().all())
    print(f"patched {len(a)} records, stock {len(b)}: identical records? {same}")


def main(run, stock_run=None):
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

    if "particle_index" in fates:
        patched_checks(fates, plan, lateral, top, end)
    if stock_run is not None:
        compare_with_stock(fates, stock_run)


if __name__ == "__main__":
    stock = sys.argv[sys.argv.index("--stock") + 1] if "--stock" in sys.argv else None
    main(sys.argv[1], stock)
