"""Compare the two 30-day runs of ``prepare_d19_comparison.py``: the effect of D19.

    scaled/20130101             every level scaled to the model surface pressure (D19)
    retrieval_levels/20130101   the retrieval's own levels

For each sounding, the column footprint is Σ_j footprint(release level j+1) × pressure weight_j ×
averaging kernel_j over the 17 released levels. Retrieval level j is matched to release level j+1
by its NUMBER, in both runs. (The legacy processing matches by nearest pressure instead,
``column.nearest_release_levels``; with D19 the releases sit at scaled pressures, so that rule
could pick a neighbouring release near the surface. New runs must match by level number.)

Printed, per sounding, scaled vs retrieval levels: ratio of column totals, of totals within 2° of
the sounding, and of the level-1 and level-2 release totals; footprint shape differences (half
L1 of the normalised maps) at native resolution and coarsened 4× and 16×. Particle noise
yardstick for 30-day runs of the same set-up (stage 1, different seeds): total scatter ~4%,
shape 32% native falling to 5% at 16×.

    uv run python dev/scripts/method_studies/compare_d19.py WORK_FOLDER
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from name_satellite.footprints import footprint_from_mixing_ratio, grid_cell_area_m2
from name_satellite.name_output import read_fields_folder
from name_satellite.validation.compare import near_field_ratio, shape_difference, total_ratio
from name_satellite.validation.legacy_pipeline import legacy_grid_coordinates

N_LEVELS = 17
RELEASE_SECONDS = 60


def read_run(folder):
    """Per-release footprints (source, lat, lon), the release plan, and the grid."""
    fields = read_fields_folder(folder / "output")
    lat = legacy_grid_coordinates(
        fields.attrs["Y grid origin"], fields.attrs["Y grid resolution"], fields.sizes["lat"]
    )
    lon = legacy_grid_coordinates(
        fields.attrs["X grid origin"], fields.attrs["X grid resolution"], fields.sizes["lon"]
    )
    footprints = footprint_from_mixing_ratio(fields.values, grid_cell_area_m2(lat, lon),
                                             RELEASE_SECONDS)  # fmt: skip
    plan = pd.read_csv(folder / "release_plan.csv", dtype={"sounding_id": str})
    if list(fields["source"].values) != plan["source_name"].tolist():
        raise ValueError(f"{folder}: fields and release plan list different sources")
    return footprints, plan, lat, lon


def per_level(footprints, plan, obs, level):
    """Footprint of release ``level`` (1-based) for each sounding of ``obs``, by sounding ID."""
    rows = plan[plan["level_number"] == level]
    position = pd.Series(rows.index, index=rows["sounding_id"])
    return footprints[position.loc[[str(i) for i in obs["sounding_id"].values]].to_numpy()]


def column(footprints, plan, obs):
    """Column footprint per sounding, retrieval level j matched to release level j+1."""
    pw = obs["pressure_weights"].values
    ak = obs["xch4_averaging_kernel"].values
    total = 0.0
    for j in range(N_LEVELS):
        weight = (pw[:, j] * ak[:, j])[:, None, None]
        total = total + per_level(footprints, plan, obs, j + 1) * weight
    return total


def describe(label, x):
    q = np.percentile(x, [5, 50, 95])
    print(f"  {label:34s} median {q[1]:.4f}   5-95% {q[0]:.4f} to {q[2]:.4f}")


def main(work):
    work = Path(work)
    scaled = work / "scaled" / "20130101"
    control = work / "retrieval_levels" / "20130101"
    with xr.open_dataset(scaled / "obs_GOSAT-BRAZIL_20130101.nc") as f:
        obs = f.load()
    fs, plan_s, lat, lon = read_run(scaled)
    fc, plan_c, lat_c, lon_c = read_run(control)
    if not (np.array_equal(lat, lat_c) and np.array_equal(lon, lon_c)):
        raise ValueError("The two runs have different grids")

    surface = pd.read_csv(scaled / "surface_pressure.csv", dtype={"sounding_id": str})
    shift = surface["model_hpa"] - surface["retrieval_hpa"]
    print(f"{len(surface)} soundings; model - retrieval surface pressure (hPa): median "
          f"{shift.median():.1f}, range {shift.min():.1f} to {shift.max():.1f}")  # fmt: skip

    col_s, col_c = column(fs, plan_s, obs), column(fc, plan_c, obs)
    print("Scaled (D19) / retrieval levels, per sounding:")
    describe("column total", total_ratio(col_s, col_c))
    near = near_field_ratio(col_s, col_c, lat, lon, obs["lat"].values, obs["lon"].values)
    describe("column total within 2°", near)
    for level in (1, 2):
        a = per_level(fs, plan_s, obs, level)
        b = per_level(fc, plan_c, obs, level)
        describe(f"level {level} release total", total_ratio(a, b))
    for factor in (1, 4, 16):
        describe(f"column shape difference, {factor}x", shape_difference(col_s, col_c, factor) / 2)

    # Does the effect follow the surface-pressure shift? (A systematic effect should.)
    ratio = total_ratio(per_level(fs, plan_s, obs, 1), per_level(fc, plan_c, obs, 1))
    by_id = surface.set_index("sounding_id").loc[[str(i) for i in obs["sounding_id"].values]]
    dp = by_id["model_hpa"] - by_id["retrieval_hpa"]
    print(f"Correlation of level-1 total ratio with the shift: "
          f"{np.corrcoef(dp.to_numpy(), ratio)[0, 1]:.2f}")  # fmt: skip


if __name__ == "__main__":
    main(sys.argv[1])
