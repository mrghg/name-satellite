"""Tests for name_satellite.release_groups (soundings sharing NAME releases, D20)."""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from name_satellite.name_input import OutputGrid
from name_satellite.release_groups import cell_size_and_origin, group_soundings, release_points

# Cells of 0.2 x 0.4 degrees with corners at whole multiples (first cell centre 0.1, 0.2).
GRID = OutputGrid(nx=10, ny=10, dx=0.4, dy=0.2, x_min=0.2, y_min=0.1)


def obs(lats, lons, times, ids=None):
    n = len(lats)
    ids = ids or [f"2013010100{i:02d}" for i in range(n)]
    return xr.Dataset(
        {
            "time": ("sounding_id", pd.to_datetime(times)),
            "lat": ("sounding_id", np.array(lats, dtype=float)),
            "lon": ("sounding_id", np.array(lons, dtype=float)),
            "pressure_levels": (
                ("sounding_id", "lev"),
                np.array([[1000.0 - i, 900.0, 800.0] for i in range(n)]),
            ),
        },
        coords={"sounding_id": ids},
    )


def test_cells_follow_the_output_grid():
    size, corner = cell_size_and_origin("grid", GRID)
    assert size == (0.2, 0.4)
    assert corner == pytest.approx((0.0, 0.0))
    assert cell_size_and_origin((1.0, 2.0), GRID)[0] == (1.0, 2.0)


def test_same_cell_and_minute_share_releases():
    o = obs(
        [0.05, 0.15, 0.25, 0.05, 0.05],  # 3rd is in the next cell north
        [0.1, 0.3, 0.1, 0.1, 0.1],
        ["2013-01-01 15:00", "2013-01-01 15:00", "2013-01-01 15:00",
         "2013-01-01 15:01", "2013-01-01 15:05"],  # 4th within 60 s, 5th a new group
    )  # fmt: skip
    assert group_soundings(o, "grid", GRID).tolist() == [0, 0, 1, 0, 2]
    assert group_soundings(o, None, GRID).tolist() == [0, 1, 2, 3, 4]
    assert group_soundings(o, (1.0, 1.0), GRID).tolist() == [0, 0, 0, 0, 1]


def test_release_point_of_a_group():
    o = obs([0.05, 0.15], [0.1, 0.3], ["2013-01-01 15:00:00", "2013-01-01 15:00:30"],
            ids=["B", "A"])  # fmt: skip
    points, members = release_points(o, np.array([0, 0]), footprint_km=10.5)
    assert points["sounding_id"].values.tolist() == ["B"]  # the earliest member's ID
    assert members["release_group"].tolist() == ["B", "B"]
    assert points["lat"].values[0] == pytest.approx(0.10)
    assert pd.Timestamp(points["time"].values[0]) == pd.Timestamp("2013-01-01 15:00:15")
    assert points["pressure_levels"].values[0].tolist() == [999.5, 900.0, 800.0]
    # Box: footprint (~0.094 deg) plus the spread of the positions.
    assert points["dlat"].values[0] == pytest.approx(0.0944 + 0.1, abs=1e-3)
    assert points["n_soundings"].values[0] == 2
    scaled = np.array([[1010.0, 900.0, 800.0], [1000.0, 900.0, 800.0]])
    points, _ = release_points(o, np.array([0, 0]), 10.5, pressure_levels_hpa=scaled)
    assert points["pressure_levels"].values[0, 0] == 1005.0
