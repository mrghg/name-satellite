"""Particle noise by release level, and the particle rates of D22 (open question 10).

Uses the stage 1 seed pair (1 Jan 2013, GOSAT-BRAZIL, 25 soundings x 17 levels, 30 days): runs B
and B2 differ only in their random seed (prepare_stage1.py, prepare_stage1_noise.py), so for every
release, (B - B2) / sqrt(2) is one run's particle noise.

For each level j, with column weight w_j = pressure weight x averaging kernel of the sounding:

* footprint: the release's footprint (mol/mol per mol/m2/s) times w_j; noise measured on the
  column total and on the cells (coarsened 4x);
* boundary: the fraction of the release's particles in each fate bin (lateral edge x 10 degrees
  along the edge x 4 height bands, the top in 10-degree latitude bands, or still inside at the
  end), times w_j. Stock 8.7.1 fate files repeat ~3% of records (X20); both runs alike.

A level's noise variance scales as 1 / (its particle count), so the noise variance of any other
allocation of particles follows, and its cost from the particles' mean lifetime (particle-hours,
from run B's fate file). Prints the shares of footprint and noise by level and compares
allocations. Caveats: one day, legacy levels (CCI v7.2, 17 levels, 19 km top); the allocations are
judged on the runs they were chosen from.

    uv run python dev/scripts/analyse_level_noise.py /work/scratch-pw5/mattrigby/name-satellite/stage1_20261109
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from name_satellite.footprints import footprint_from_mixing_ratio, grid_cell_area_m2
from name_satellite.name_output import read_fields, read_name_fate_file
from name_satellite.refdata import ref_path
from name_satellite.validation import coarsen

LEVELS = 17  # released by the stage 1 runs
COMPARE = 16  # levels modelled by new runs (D21)
SOUNDINGS = 25
LEGACY_RATES = np.array([30] * 12 + [3] * 5)  # particles/s
RELEASE_SECONDS = 60

# Stage 1 dispersion domain (input.txt) and binning of the particle fates.
LON_MIN, LON_MAX, LAT_MIN, LAT_MAX, TOP_M = -96.329, -19.449, -64.981, 27.557, 19000.0
ALONG_EDGE_DEG = 10.0
HEIGHT_EDGES_M = [2000.0, 5000.0, 10000.0]
MAX_AGE_H = 720.0

ALLOCATIONS = {  # particles/s for levels 1..16
    "legacy (30/s to 12, 3/s above)": [30] * 12 + [3] * 4,
    "13 at 30/s": [30] * 13 + [3] * 3,
    "13 at 30/s, 14 at 10/s (D22)": [30] * 13 + [10, 3, 3],
    "13-14 at 30/s": [30] * 14 + [3] * 2,
    "all at 30/s": [30] * 16,
}


def column_weights():
    """(sounding, level) pressure weight x averaging kernel from the legacy observation file."""
    obs = xr.open_dataset(ref_path("obs"))
    return (obs["pressure_weights"].values * obs["xch4_averaging_kernel"].values)[:, :LEVELS]


def weighted_footprints(run, weights):
    """(sounding, level, lat, lon) footprints of one run, times the column weights."""
    fields = read_fields(sorted((run / "output").glob("Fields_grid1_C1_T1_*.txt"))[0])
    area = grid_cell_area_m2(fields["lat"].values, fields["lon"].values)
    fp = footprint_from_mixing_ratio(fields.values, area, float(RELEASE_SECONDS))
    fp = fp.reshape(SOUNDINGS, LEVELS, *area.shape)  # sources are sounding-major, levels within
    return fp * weights[:, :, None, None]


def fate_bin_counts(run):
    """Particles per (source, fate bin) of one run, and the particle table."""
    particles = read_name_fate_file(run / "output" / "Particles_C1.txt")[0]
    distance = pd.DataFrame(
        {
            "W": particles["lon"] - LON_MIN,
            "E": LON_MAX - particles["lon"],
            "S": particles["lat"] - LAT_MIN,
            "N": LAT_MAX - particles["lat"],
        }
    )
    edge = distance.idxmin(axis=1)  # nearest lateral face
    along = np.where(edge.isin(["W", "E"]), particles["lat"] - LAT_MIN, particles["lon"] - LON_MIN)
    along_bin = (along // ALONG_EDGE_DEG).astype(int).astype(str)
    height_bin = np.digitize(particles["z"], HEIGHT_EDGES_M).astype(str)
    top_bin = ((particles["lat"] - LAT_MIN) // ALONG_EDGE_DEG).astype(int).astype(str)
    lateral = edge + "_" + along_bin + "_" + height_bin
    exit_bin = np.where(particles["z"] >= TOP_M - 50.0, "top_" + top_bin, lateral)
    fate = np.where(particles["age_hours"] >= MAX_AGE_H - 0.1, "inside", exit_bin)
    return pd.crosstab(particles["source_number"], fate), particles


def weighted_fate_fractions(counts, columns, weights):
    """(sounding, level, bin) fractions of the particles released, times the column weights."""
    counts = counts.reindex(columns=columns, fill_value=0)
    counts = counts.reindex(range(1, SOUNDINGS * LEVELS + 1), fill_value=0)
    released = np.tile(LEGACY_RATES * RELEASE_SECONDS, SOUNDINGS)[:, None]
    fractions = (counts.to_numpy() / released).reshape(SOUNDINGS, LEVELS, -1)
    return fractions * weights[:, :, None]


def main(work):
    work = Path(work)
    weights = column_weights()

    fb = weighted_footprints(work / "run_B", weights)
    f2 = weighted_footprints(work / "run_B2", weights)
    mean_fp = (fb + f2) / 2
    noise_fp = (fb - f2) / np.sqrt(2)
    total_var = (noise_fp.sum(axis=(2, 3)) ** 2).mean(axis=0)
    cells = np.stack([coarsen(noise_fp[:, j], 4) for j in range(LEVELS)], axis=1)
    cell_var = (cells**2).sum(axis=(2, 3)).mean(axis=0)

    counts_b, particles_b = fate_bin_counts(work / "run_B")
    counts_2, _ = fate_bin_counts(work / "run_B2")
    columns = sorted(set(counts_b.columns) | set(counts_2.columns))
    hb = weighted_fate_fractions(counts_b, columns, weights)
    h2 = weighted_fate_fractions(counts_2, columns, weights)
    boundary_var = ((((hb - h2) / np.sqrt(2)) ** 2).sum(axis=2)).mean(axis=0)

    level = (particles_b["source_number"] - 1) % LEVELS + 1
    mean_age = particles_b.groupby(level)["age_hours"].mean().to_numpy()

    k = slice(0, COMPARE)
    share = mean_fp.sum(axis=(0, 2, 3))[k] / mean_fp.sum()
    table = pd.DataFrame(
        {
            "particles/s": LEGACY_RATES[k],
            "weight share": np.median(weights[:, k] / weights[:, k].sum(1, keepdims=True), 0),
            "footprint share": share,
            "noise: column total": total_var[k] / total_var[k].sum(),
            "noise: cells 4x": cell_var[k] / cell_var[k].sum(),
            "noise: boundary": boundary_var[k] / boundary_var[k].sum(),
            "mean age (h)": mean_age[k],
        },
        index=pd.RangeIndex(1, COMPARE + 1, name="level"),
    )
    print(f"{len(columns)} fate bins. Shares by level (stage 1 seed pair):")
    print(table.round(3).to_string())

    # Noise variance of level j scales as 1 / particles; cost as particles x mean age.
    rows = []
    legacy = LEGACY_RATES[k].astype(float)
    cost0 = (legacy * mean_age[k]).sum()
    for name, rates in ALLOCATIONS.items():
        rates = np.asarray(rates, dtype=float)
        row = {"allocation": name, "cost": (rates * mean_age[k]).sum() / cost0 - 1}
        for label, var in [("column total", total_var), ("cells 4x", cell_var),
                           ("boundary", boundary_var)]:  # fmt: skip
            row[label] = (var[k] * legacy / rates).sum() / var[k].sum()
        rows.append(row)
    print("\nNoise variance relative to legacy, and extra cost (particle-hours):")
    print(pd.DataFrame(rows).set_index("allocation").round(2).to_string())


if __name__ == "__main__":
    main(sys.argv[1])
