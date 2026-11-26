"""Reading the future-method month file (future.footprint_file), on a small synthetic file."""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from name_satellite.footprints import grid_cell_area_m2
from name_satellite.future.footprint_file import (
    REQUIRED_VARIABLES,
    aggregate,
    box_edges,
    finest_footprint,
    footprint_density,
    open_footprints,
    sounding_footprint,
    sounding_table,
)
from name_satellite.future.grids import (
    FINE,
    MEDIUM,
    box_offset,
    box_shape,
    global_centres,
    snapped_box,
)

SOUNDINGS = {"s1": (-60.3, -10.2), "s2": (-47.9, -23.6)}  # sounding ID: (lon, lat)


def block_sum(values, factor):
    """Sum factor x factor blocks of a (ny, nx) array."""
    ny, nx = values.shape
    return values.reshape(ny // factor, factor, nx // factor, factor).sum(axis=(1, 3))


def place_in_block(values, n):
    out = np.zeros((n, n))
    out[: values.shape[0], : values.shape[1]] = values
    return out


@pytest.fixture(scope="module")
def synthetic():
    """A file with two soundings whose three levels come from one 0.1° global field each, so
    the levels agree exactly. Returns (dataset, {sounding: 0.1° global field})."""
    rng = np.random.default_rng(1)
    truth, parts = {}, {"fine": [], "medium": [], "coarse": []}
    offsets, shapes = {"fine": [], "medium": []}, {"fine": [], "medium": []}
    for sounding, (lon, lat) in SOUNDINGS.items():
        field = rng.random((1800, 3600))
        truth[sounding] = field
        parts["coarse"].append(block_sum(field, 25))
        for level, factor in ((FINE, 1), (MEDIUM, 5)):
            box = snapped_box(lon, lat, level.half_width_deg)
            (j0, i0), (ny, nx) = box_offset(box, level.cell_deg), box_shape(box, level.cell_deg)
            on_level = block_sum(field, factor)
            parts[level.name].append(place_in_block(on_level[j0 : j0 + ny, i0 : i0 + nx],
                                                    level.block_cells))  # fmt: skip
            offsets[level.name].append((j0, i0))
            shapes[level.name].append((ny, nx))
    fine_lat, fine_lon = global_centres(0.1)
    medium_lat, medium_lon = global_centres(0.5)
    coarse_lat, coarse_lon = global_centres(2.5)
    s = "sounding_id"
    ids = list(SOUNDINGS)
    times = pd.to_datetime(["2013-01-01 14:00", "2013-01-01 15:00"]).values
    ds = xr.Dataset(
        {
            "srr_fine": ((s, "fine_y", "fine_x"), np.array(parts["fine"])),
            "fine_offset": ((s, "yx"), np.array(offsets["fine"])),
            "fine_shape": ((s, "yx"), np.array(shapes["fine"])),
            "srr_medium": ((s, "medium_y", "medium_x"), np.array(parts["medium"])),
            "medium_offset": ((s, "yx"), np.array(offsets["medium"])),
            "medium_shape": ((s, "yx"), np.array(shapes["medium"])),
            "srr_coarse": ((s, "coarse_lat", "coarse_lon"), np.array(parts["coarse"])),
            "end_positions": (
                (s, "end_height", "coarse_lat", "coarse_lon"),
                np.full((2, 21, 72, 144), 0.5 / (21 * 72 * 144)),
            ),
            "column_weight": ((s, "release_level"), np.full((2, 16), 1 / 32)),
            "time": (s, times),
            "lat": (s, [v[1] for v in SOUNDINGS.values()]),
            "lon": (s, [v[0] for v in SOUNDINGS.values()]),
            "release_group": (s, np.array([0, 1])),
            "release_lat": (s, [v[1] for v in SOUNDINGS.values()]),
            "release_lon": (s, [v[0] for v in SOUNDINGS.values()]),
            "release_time": (s, times),
            "end_time": (s, times - np.timedelta64(240, "h")),
            "top_exit_mean_age": (s, [100.0, 120.0]),
        },
        coords={
            s: ids,
            "release_level": np.arange(1, 17),
            "fine_lat_global": fine_lat,
            "fine_lon_global": fine_lon,
            "medium_lat_global": medium_lat,
            "medium_lon_global": medium_lon,
            "coarse_lat": coarse_lat,
            "coarse_lon": coarse_lon,
            "end_height": np.append(np.arange(500.0, 20000.0, 1000.0), np.nan),
        },
    )
    return ds, truth


def test_open_refuses_a_file_without_the_layout(tmp_path, synthetic):
    ds, _ = synthetic
    good = tmp_path / "good.nc"
    ds.isel(sounding_id=[0]).to_netcdf(good)
    assert set(REQUIRED_VARIABLES) <= set(open_footprints(good).variables)
    bad = tmp_path / "bad.nc"
    ds.drop_vars("end_time").isel(sounding_id=[0]).to_netcdf(bad)
    with pytest.raises(ValueError, match="end_time"):
        open_footprints(bad)


def test_sounding_footprint_is_the_box_cut_from_the_global_grid(synthetic):
    ds, truth = synthetic
    fine = sounding_footprint(ds, "s1", "fine")
    box = box_edges(ds, "s1", "fine")
    assert box == snapped_box(*SOUNDINGS["s1"], 5.0)
    j0, i0 = box_offset(box, 0.1)
    ny, nx = box_shape(box, 0.1)
    np.testing.assert_allclose(fine.values, truth["s1"][j0 : j0 + ny, i0 : i0 + nx])
    assert fine["lat"].values[0] == pytest.approx(box.south + 0.05)
    assert fine["lon"].values[-1] == pytest.approx(box.east - 0.05)
    # 0.1° cells summed to 0.5° are the 0.5° cells there
    medium = sounding_footprint(ds, "s1", "medium")
    summed = aggregate(fine, 5)
    np.testing.assert_allclose(summed.values, medium.sel(lat=summed["lat"], lon=summed["lon"],
                                                         method="nearest").values)  # fmt: skip


def test_finest_footprint_takes_the_finest_level_and_spreads_the_coarser(synthetic):
    ds, truth = synthetic
    one = finest_footprint(ds, "fine", "s1").values
    field = truth["s1"]
    fine_box, medium_box = snapped_box(*SOUNDINGS["s1"], 5.0), snapped_box(*SOUNDINGS["s1"], 25.0)
    (j0, i0), (ny, nx) = box_offset(fine_box, 0.1), box_shape(fine_box, 0.1)
    np.testing.assert_allclose(one[j0 : j0 + ny, i0 : i0 + nx], field[j0 : j0 + ny, i0 : i0 + nx])
    # Inside the medium box, outside the fine box: each 0.5° cell spread over its 25 cells
    mj0, mi0 = box_offset(medium_box, 0.1)
    np.testing.assert_allclose(one[mj0 : mj0 + 5, mi0 : mi0 + 5],
                               field[mj0 : mj0 + 5, mi0 : mi0 + 5].sum() / 25)  # fmt: skip
    # Outside the medium box: each 2.5° cell spread over its 625 cells
    np.testing.assert_allclose(one[:25, :25], field[:25, :25].sum() / 625)
    np.testing.assert_allclose(one.sum(), field.sum())
    # A sum of soundings is the sum of their footprints, on every level
    for level in ("coarse", "medium", "fine"):
        both = finest_footprint(ds, level)
        each = finest_footprint(ds, level, "s1") + finest_footprint(ds, level, "s2")
        np.testing.assert_allclose(both.values, each.values)


def test_sounding_table_fractions_in_the_boxes(synthetic):
    ds, truth = synthetic
    table = sounding_table(ds)
    assert list(table.index) == ["s1", "s2"]
    for sounding, (lon, lat) in SOUNDINGS.items():
        field = truth[sounding]
        for name, half in (("fine", 5.0), ("medium", 25.0)):
            box = snapped_box(lon, lat, half)
            (j0, i0), (ny, nx) = box_offset(box, 0.1), box_shape(box, 0.1)
            expected = field[j0 : j0 + ny, i0 : i0 + nx].sum() / field.sum()
            assert table.loc[sounding, f"in_{name}_box"] == pytest.approx(expected)
    assert table["column_weight"].tolist() == pytest.approx([0.5, 0.5])
    assert table["left_top"].tolist() == pytest.approx([1 / 21, 1 / 21])


def test_density_is_the_same_on_every_level_for_a_uniform_density(synthetic):
    ds, _ = synthetic
    lat, lon = global_centres(0.1)
    area = grid_cell_area_m2(lat, lon)
    uniform = ds.copy()
    uniform["srr_coarse"] = (("sounding_id", "coarse_lat", "coarse_lon"),
                             np.broadcast_to(block_sum(area, 25), (2, 72, 144)))  # fmt: skip
    j0, i0 = ds["medium_offset"].values[0]
    ny, nx = ds["medium_shape"].values[0]
    medium = place_in_block(block_sum(area, 5)[j0 : j0 + ny, i0 : i0 + nx], 105)
    uniform["srr_medium"] = (("sounding_id", "medium_y", "medium_x"), np.array([medium, medium]))
    for level in ("coarse", "medium"):
        density = footprint_density(sounding_footprint(uniform, "s1", level))
        np.testing.assert_allclose(density.values, 1.0, rtol=1e-3)
