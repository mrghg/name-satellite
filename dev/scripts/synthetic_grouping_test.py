"""Synthetic test of D20: can a group's release stand in for each member's own release?

D20 lets soundings in the same output-grid cell and minute share one set of NAME releases, at the
group's mean position, in a box of the largest member footprint plus the members' spread. A
member's own release would differ in two ways, both far below the grid spacing:

1. its centre is offset from the group's by up to ~2 km (January 2013);
2. its box is narrower (its own ~10.5 km footprint, without the spread).

Instead of rerunning NAME, this script applies those two differences to the processed footprints
of the grouped run, assuming the flow is locally uniform over a few km (so a small change of the
release becomes the same change of the footprint):

* offset: shift the footprint by the member's offset, by linear interpolation between neighbouring
  cells (a sub-cell shift);
* box: a release spread uniformly over a box of width w adds variance w²/12 along that axis, so
  narrowing the box from the group's to the member's removes variance (w_group² − w_own²)/12.
  Applied (with the opposite sign the effect is the same size) as a three-point blur with that
  variance, in grid cells².

The difference made is then compared with pure particle noise: the stage 1 seed pair (two runs of
1 Jan 2013 that differ only in their random seed), with the same measures
(``validation.shape_difference`` at several coarsenings, ``validation.total_ratio``). The seed
pair used the legacy particle rates, so its noise is somewhat larger than the new runs' (D22:
variance −37% in cells); both are printed.

Usage:
    uv run python dev/scripts/synthetic_grouping_test.py PROCESSED_FILE OBS_FILE STAGE1_FOLDER
"""

import sys

import numpy as np
import xarray as xr

from name_satellite.releases import footprint_box_degrees
from name_satellite.validation import shape_difference, total_ratio

FOOTPRINT_KM = 10.5  # GOSAT-1 (domains.GOSAT_BRAZIL.footprint_km)
COARSENINGS = (1, 2, 4)


def shift_axis(f, cells, axis):
    """Shift 2-D field f by a fraction of a cell (|cells| <= 1) along axis, linearly.

    Content moves towards higher index for cells > 0. Cells shifted in from outside are 0.
    """
    if abs(cells) > 1:
        raise ValueError(f"Shift of {cells} cells; only sub-cell shifts are supported")
    step = 1 if cells > 0 else -1
    neighbour = np.roll(f, step, axis=axis)
    edge = [slice(None)] * 2
    edge[axis] = 0 if step == 1 else -1
    neighbour[tuple(edge)] = 0.0  # don't wrap around the grid
    return (1 - abs(cells)) * f + abs(cells) * neighbour


def blur_axis(f, variance_cells2, axis):
    """Three-point blur along axis with the given variance (cells², at most 0.5)."""
    if not 0 <= variance_cells2 <= 0.5:
        raise ValueError(f"Blur variance {variance_cells2} outside [0, 0.5] cells²")
    up = np.roll(f, 1, axis=axis)
    down = np.roll(f, -1, axis=axis)
    return (1 - variance_cells2) * f + variance_cells2 / 2 * (up + down)


def as_member(fp, dlat_cells, dlon_cells, var_lat_cells2, var_lon_cells2):
    """The group footprint fp (lat, lon) changed as if released from the member's own box."""
    out = shift_axis(fp, dlat_cells, axis=0) if dlat_cells else fp
    out = shift_axis(out, dlon_cells, axis=1) if dlon_cells else out
    out = blur_axis(out, var_lat_cells2, axis=0)
    return blur_axis(out, var_lon_cells2, axis=1)


def main(processed_file, obs_file, stage1_folder):
    ds = xr.open_dataset(processed_file).load()
    obs = xr.open_dataset(obs_file).load()
    lat, lon = ds["latitude"].values, ds["longitude"].values
    cell_lat, cell_lon = np.diff(lat).mean(), np.diff(lon).mean()
    srr = ds["srr"].values[ds["footprint_index"].values]  # (sounding, lat, lon)

    # Member positions, by sounding ID.
    obs_lat = dict(zip(obs["sounding_id"].values.astype(str), obs["lat"].values, strict=True))
    obs_lon = dict(zip(obs["sounding_id"].values.astype(str), obs["lon"].values, strict=True))
    groups = ds["release_group"].values
    group_size = {g: (groups == g).sum() for g in np.unique(groups)}

    own, moved, offsets_km = [], [], []
    for i, sid in enumerate(ds["sounding_id"].values.astype(str)):
        if group_size[groups[i]] < 2:
            continue
        dlat = obs_lat[sid] - float(ds["release_lat"][i])
        dlon = obs_lon[sid] - float(ds["release_lon"][i])
        own_dlon, own_dlat = footprint_box_degrees(obs_lat[sid], FOOTPRINT_KM)
        var_lat = (float(ds["release_dlat"][i]) ** 2 - own_dlat**2) / 12 / cell_lat**2
        var_lon = (float(ds["release_dlon"][i]) ** 2 - own_dlon**2) / 12 / cell_lon**2
        own.append(srr[i])
        moved.append(as_member(srr[i], dlat / cell_lat, dlon / cell_lon, var_lat, var_lon))
        offsets_km.append(np.hypot(dlat * 111.2, dlon * 111.2 * np.cos(np.radians(obs_lat[sid]))))
    own, moved, offsets_km = np.array(own), np.array(moved), np.array(offsets_km)
    print(f"{len(own)} soundings in multi-member groups; offset from the group's release centre "
          f"median {np.median(offsets_km):.2f} km, max {offsets_km.max():.2f} km "
          f"(grid cell {cell_lat * 111.2:.0f} x {cell_lon * 111.2:.0f} km at the equator)")  # fmt: skip

    # Particle noise: stage 1 seed pair (fp is (lat, lon, time) in the legacy layout).
    b = xr.open_dataset(f"{stage1_folder}/run_B/fp_20130101_legacy_processing.nc")
    b2 = xr.open_dataset(f"{stage1_folder}/run_B2/fp_20130101_legacy_processing.nc")
    fb = b["fp"].transpose("time", "lat", "lon").values.astype(float)
    fb2 = b2["fp"].transpose("time", "lat", "lon").values.astype(float)

    def summary(values):
        return f"median {np.median(values):.2e}  90% {np.percentile(values, 90):.2e}  max {values.max():.2e}"

    print("\nTotal, |ratio − 1|:")
    print("  member release vs group release:", summary(np.abs(total_ratio(moved, own) - 1)))
    print("  seed pair (particle noise):     ", summary(np.abs(total_ratio(fb2, fb) - 1)))
    print("\nShape difference (L1 of normalised footprints, 0..2):")
    for factor in COARSENINGS:
        effect = shape_difference(moved, own, factor)
        noise = shape_difference(fb2, fb, factor)
        print(f"  {factor}x{factor} cells: member vs group {summary(effect)}")
        print(f"             seed pair      {summary(noise)}"
              f"   (D22-scaled median {np.median(noise) * np.sqrt(1 - 0.37):.2e})")  # fmt: skip
        print(
            f"             ratio of medians, effect / noise: {np.median(effect) / np.median(noise):.3f}"
        )


if __name__ == "__main__":
    main(*sys.argv[1:4])
