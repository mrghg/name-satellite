"""GLIDE spec §10 checks that apply to whole delivered blocks (not only the one-day check).

For each delivered footprint file (block or single-release run) and its ``_end`` file in FOLDER, prints: releases, particles
per release (§10 check 7), whether there is one end record per particle, the oldest particle
(check 2: at most 240 h plus one 300 s step), footprint totals (check 3, m2 s/mol) and the
fractions of particles that left through the sides and top.

    uv run python dev/scripts/glide_tests/check_blocks.py FOLDER
"""

import sys
from pathlib import Path

import numpy as np
import xarray as xr

from name_satellite.glide_runs import spec

MAX_AGE = spec.MAX_AGE_HOURS + spec.SYNC_MINUTES / 60.0  # h


def check(footprint_path):
    ds = xr.open_dataset(footprint_path)
    end = xr.open_dataset(str(footprint_path).replace(".nc", "_end.nc"))
    n = ds["n_particles"].values
    rows = end["time"].to_series().value_counts().reindex(ds["time"].values).fillna(0).values
    totals = ds["srr"].sum(["latitude", "longitude"]).values
    oldest = float(end["age"].max())
    ok = np.array_equal(rows, n) and oldest <= MAX_AGE
    print(
        f"{'PASS' if ok else 'FAIL'} {footprint_path.name}: {ds.sizes['time']} releases; "
        f"n_particles {n.min()}-{n.max()}; end records match: {np.array_equal(rows, n)}; "
        f"oldest {oldest:.4f} h; totals min {totals.min():.2f} median {np.median(totals):.2f} "
        f"max {totals.max():.2f} m2 s/mol ({((totals >= 5) & (totals <= 17)).mean():.0%} in 5-17); "
        f"left sides {ds['fraction_left_sides'].values.mean():.2f}, "
        f"top {ds['fraction_left_top'].values.mean():.3f} (block means)"
    )


if __name__ == "__main__":
    for path in sorted(Path(sys.argv[1]).glob("*.nc")):
        if not path.name.endswith(("_end.nc", "_positions.nc")):  # footprint files only
            check(path)
