"""Compare a newly processed footprint file with a legacy reference file.

Used after dev/scripts/legacy_validation/run_legacy_process.py: does today's legacy processing code reproduce the
footprint files that were made in 2018/2019?

Only the soundings present in the NEW file are compared (it may cover just a few days).
Soundings are matched by time; the files are read-only here.

    uv run python dev/scripts/legacy_validation/compare_processed.py NEW_FILE REFERENCE_FILE

For each variable it reports whether values are identical, and otherwise the largest
relative difference. The footprint is called `fp` in old files and `srr` in newer ones; both
names are tried.
"""

import sys

import numpy as np
import xarray as xr

# Old (2018/2019 "fp" format) name  ->  newer ("srr" format) name.
NAME_PAIRS = {
    "fp": "srr",
    "lat": "latitude",
    "lon": "longitude",
    "temperature": "air_temperature",
    "pressure": "air_pressure",
    "wind_speed": "wind_speed",
    "wind_direction": "wind_from_direction",
    "PBLH": "atmosphere_boundary_layer_thickness",
    "release_lon": "release_lon",
    "release_lat": "release_lat",
    "particle_locations_n": "particle_locations_n",
    "particle_locations_e": "particle_locations_e",
    "particle_locations_s": "particle_locations_s",
    "particle_locations_w": "particle_locations_w",
}


def get(ds: xr.Dataset, old_name: str) -> xr.DataArray:
    """Fetch a variable by its old name, or its newer equivalent; rename dims to old names."""
    for name in (old_name, NAME_PAIRS[old_name]):
        if name in ds.variables:
            da = ds[name]
            return da.rename({d: "lat" for d in da.dims if d == "latitude"}).rename(
                {d: "lon" for d in da.dims if d == "longitude"}
            )
    raise KeyError(f"{old_name} / {NAME_PAIRS[old_name]} not in file")


def describe_difference(new: np.ndarray, ref: np.ndarray) -> str:
    """One-line summary of how two arrays of the same shape differ."""
    if new.shape != ref.shape:
        return f"SHAPES DIFFER new {new.shape} vs reference {ref.shape}"
    if np.array_equal(new, ref, equal_nan=True):
        return "identical"
    scale = np.nanmax(np.abs(ref))
    max_abs = np.nanmax(np.abs(new - ref))
    text = f"max |diff| = {max_abs:.3g} ({max_abs / scale:.2e} of max |reference|)"
    # Hint at pure unit changes: a constant offset (e.g. C -> K) or constant factor (Pa -> hPa).
    offset = new - ref
    if np.nanmax(offset) - np.nanmin(offset) < 1e-3 * max(scale, 1):
        text += f"; constant offset {np.nanmean(offset):.4g}"
    with np.errstate(divide="ignore", invalid="ignore"):
        factor = np.where(ref != 0, new / ref, np.nan)
    if np.nanmax(factor) - np.nanmin(factor) < 1e-4 * abs(np.nanmean(factor)):
        text += f"; constant factor {np.nanmean(factor):.6g}"
    return text


def main(new_path: str, ref_path: str) -> None:
    new = xr.open_dataset(new_path)
    ref = xr.open_dataset(ref_path)

    # Match soundings by time.
    missing = np.setdiff1d(new.time.values, ref.time.values)
    if missing.size:
        raise ValueError(f"{missing.size} times in the new file are not in the reference")
    ref = ref.sel(time=new.time.values)
    print(f"Comparing {new.sizes['time']} soundings")
    for coord in ("lat", "lon"):
        print(f"{coord:22s}", describe_difference(get(new, coord).values, get(ref, coord).values))

    # Footprint: put both in the same orientation before comparing.
    fp_new = get(new, "fp").transpose("lat", "lon", "time").values
    fp_ref = get(ref, "fp").transpose("lat", "lon", "time").values
    print(f"{'footprint':22s}", describe_difference(fp_new, fp_ref))
    ratio = fp_new.sum(axis=(0, 1)) / fp_ref.sum(axis=(0, 1))
    print(
        f"{'':22s} total ratio new/ref per sounding: min {ratio.min():.6f}, max {ratio.max():.6f}"
    )

    for name in list(NAME_PAIRS)[3:]:  # skip fp, lat, lon (done above)
        try:
            a, b = get(new, name), get(ref, name)
        except KeyError as err:
            print(f"{name:22s} {err}")
            continue
        print(f"{name:22s}", describe_difference(a.transpose(*b.dims).values, b.values))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
