"""Tests for name_satellite.validation (footprint comparison measures)."""

import numpy as np

from name_satellite.validation import (
    coarsen,
    histogram_difference,
    near_field_ratio,
    shape_difference,
    total_ratio,
)


def test_total_and_shape():
    g = np.zeros((1, 4, 4))
    g[0, 0, 0] = 1.0
    f = 3 * g
    assert total_ratio(f, g).tolist() == [3.0]
    assert shape_difference(f, g).tolist() == [0.0]  # same shape, different total
    h = np.zeros((1, 4, 4))
    h[0, 3, 3] = 1.0
    assert shape_difference(h, g).tolist() == [2.0]  # no overlap
    assert shape_difference(h, g, factor=4).tolist() == [0.0]  # same 4x4 block


def test_coarsen_drops_leftover_cells():
    f = np.ones((2, 5, 7))
    c = coarsen(f, 2)
    assert c.shape == (2, 2, 3)
    assert (c == 4).all()


def test_near_field_ratio():
    lat = np.array([0.0, 1.0, 5.0])
    lon = np.array([0.0, 1.0, 5.0])
    g = np.ones((1, 3, 3))
    f = g.copy()
    f[0, 2, 2] = 10.0  # far from the sounding at (0, 0): not counted
    f[0, 1, 1] = 2.0
    assert near_field_ratio(f, g, lat, lon, [0.0], [0.0]).tolist() == [5.0 / 4.0]


def test_histogram_difference():
    a = {e: np.ones((1, 2, 3)) for e in "nesw"}
    b = {e: np.ones((1, 2, 3)) for e in "nesw"}
    assert histogram_difference(a, b).tolist() == [0.0]
    a["n"] = a["n"] * 2
    assert histogram_difference(a, b).tolist() == [6.0 / 24.0]
