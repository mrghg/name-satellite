"""Compare name_satellite.observations with the legacy GOSAT preprocessing on the same L2 file.

The legacy reference is written by run_legacy_gosat.py (binned.nc and unbinned.nc). Soundings are
matched by exposure ID (for bins: the comma-joined IDs of the bin's soundings).

    uv run python dev/scripts/legacy_validation/compare_observations_with_legacy.py L2_FILE LEGACY_FOLDER
"""

import sys

import numpy as np
import xarray as xr

from name_satellite.observations import (
    PROFILE_VARIABLES,
    bad_pressure_mask,
    latlon_mask,
    legacy_bin,
    quality_mask,
    read_l2,
    select,
)

LAT_BOUNDS = (-35.753, 7.251)  # as run_legacy_gosat.py
LON_BOUNDS = (-75.984, -32.793)
BIN_SIZE = (0.234, 0.352)
COMPARED = [
    "latitude",
    "longitude",
    "xch4",
    "xch4_uncertainty",
    "pressure_levels",
    "pressure_weight",
    "xch4_averaging_kernel",
    "ch4_profile_apriori",
]


def compare(ours, legacy, label, extra=()):
    legacy_ids = list(legacy["exposure_id"].values)
    our_ids = list(ours["exposure_id"].values)
    print(f"{label}: {len(our_ids)} ours, {len(legacy_ids)} legacy, "
          f"same ID sets: {set(our_ids) == set(legacy_ids)}")  # fmt: skip
    if set(our_ids) != set(legacy_ids):
        return
    order = [our_ids.index(i) for i in legacy_ids]
    ours = ours.isel(sounding=order)
    for name in [*COMPARED, *extra]:
        a = ours[name].values.astype(float)
        b = legacy[name].values.astype(float)
        relative = np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-30))
        print(f"   {name:24s} max relative difference {relative:.1e}")
    dt = (ours["time"].values - legacy["time"].values) / np.timedelta64(1, "us")
    print(f"   {'time':24s} max |difference| {np.abs(dt).max():.0f} µs")


def main(l2_file, legacy_folder):
    table = read_l2(l2_file)
    table = select(table, quality_mask(table) & bad_pressure_mask(table))
    print(f"{table.attrs['product']} {table.attrs['product_version']}; "
          f"time resolution: {table.attrs['time_resolution']}")  # fmt: skip

    unbinned = select(table, latlon_mask(table, LAT_BOUNDS, LON_BOUNDS))
    compare(unbinned, xr.open_dataset(f"{legacy_folder}/unbinned.nc"), "unbinned")

    binned = legacy_bin(table, LAT_BOUNDS, LON_BOUNDS, BIN_SIZE)
    compare(binned, xr.open_dataset(f"{legacy_folder}/binned.nc"), "binned", extra=("dlat", "dlon"))
    assert set(PROFILE_VARIABLES) <= set(COMPARED)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
