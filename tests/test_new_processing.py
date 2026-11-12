"""Tests of the processing of new runs (processing.py): particle fates, columns, the file."""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from name_satellite.name_input import DispersionDomain
from name_satellite.processing import (
    BITROUND_BITS,
    HEIGHT_BIN_EDGES_M,
    _nearest_cell,
    assemble,
    classify_fates,
    release_fates,
    write_footprint_file,
)

DOMAIN = DispersionDomain(x_min=-10.0, x_max=10.0, y_min=-5.0, y_max=5.0, z_max_m_agl=20000.0)
LAT = np.array([-4.0, -2.0, 0.0, 2.0, 4.0])
LON = np.array([-8.0, -4.0, 0.0, 4.0, 8.0])


def fate_table(rows):
    """Fate records from (source, index, age_hours, lon, lat, z, code) tuples."""
    columns = ["source_number", "particle_index", "age_hours", "lon", "lat", "z", "fate_code"]
    return pd.DataFrame(rows, columns=columns)


def test_classify_fates():
    particles = fate_table(
        [
            (1, 1, 10.0, 0.0, 5.01, 500.0, 8),  # north
            (1, 2, 11.0, 10.2, 0.0, 1500.0, 8),  # east
            (1, 3, 12.0, 0.0, -5.1, 2500.0, 8),  # south
            (1, 4, 13.0, -10.05, 0.0, 3500.0, 8),  # west
            (1, 5, 14.0, 0.0, 0.0, 20012.0, 8),  # top
            (1, 6, 15.0, 10.3, 5.1, 900.0, 8),  # corner: further beyond east (0.3 > 0.1)
            (1, 7, 16.0, 10.3, 0.0, 20005.0, 8),  # beyond a side and above the top: top
            (1, 8, 720.0, 1.0, 1.0, 800.0, 16),  # still inside
            (1, 9, 720.0, 10.02, 1.0, 200.0, 24),  # left on its last step (8 + 16): east
            (1, 10, 50.0, 10.0, 1.0, 700.0, 8),  # on the east edge as printed: east
            (1, 11, 50.0, 1.0, 1.0, 20000.0, 8),  # at the top as printed: top
        ]
    )
    assert classify_fates(particles, DOMAIN).tolist() == [
        "n", "e", "s", "w", "top", "e", "top", "inside", "e", "e", "top"
    ]  # fmt: skip


def test_classify_fates_rejects_unexplained_records():
    for code in (8, 24):
        inside_but_left = fate_table([(1, 1, 3.0, 0.0, 0.0, 500.0, code)])
        with pytest.raises(ValueError, match="inside it"):
            classify_fates(inside_but_left, DOMAIN)
    outside_but_inside = fate_table([(1, 1, 720.0, 11.0, 0.0, 500.0, 16)])
    with pytest.raises(ValueError, match="outside the domain"):
        classify_fates(outside_but_inside, DOMAIN)
    unknown_code = fate_table([(1, 1, 3.0, 11.0, 0.0, 500.0, 4)])
    with pytest.raises(ValueError, match="fate codes"):
        classify_fates(unknown_code, DOMAIN)


def test_nearest_cell_clips_to_the_grid():
    assert _nearest_cell(np.array([-9.9, -6.1, 1.9, 2.1, 30.0]), LON).tolist() == [0, 0, 2, 3, 4]


def test_release_fates_fractions_and_ages():
    particles = fate_table(
        [
            (1, 1, 10.0, 0.0, 5.01, 500.0, 8),  # north, height bin 0, lon cell 2
            (1, 2, 30.0, 0.1, 5.02, 999.0, 8),  # same bin
            (1, 3, 12.0, 10.2, -3.9, 20000.0, 8),  # east, top height bin, lat cell 0
            (1, 4, 14.0, 7.9, 2.1, 20012.0, 8),  # top, cell (lat 3, lon 4)
            (2, 1, 720.0, 0.0, 0.0, 100.0, 16),  # source 2: inside
            (2, 2, 5.0, -10.1, 0.0, 1000.0, 8),  # west, height bin 1
        ]
    )
    fate = classify_fates(particles, DOMAIN)
    result = release_fates(particles, fate, [1, 2], LAT, LON)
    assert result["n_released"].tolist() == [4, 2]
    assert result["n"].shape == (2, len(HEIGHT_BIN_EDGES_M) - 1, len(LON))
    assert result["n"][0, 0, 2] == 0.5
    assert result["age_n"][0, 0, 2] / result["n"][0, 0, 2] == 20.0  # mean of 10 h and 30 h
    assert result["e"][0, 19, 0] == 0.25
    assert result["top"][0, 3, 4] == 0.25
    assert result["inside"].tolist() == [0.0, 0.5]
    assert result["w"][1, 1, 2] == 0.5
    totals = sum(result[f].sum(axis=(1, 2)) for f in ("n", "e", "s", "w", "top"))
    np.testing.assert_allclose(totals + result["inside"], 1.0)


def test_release_fates_needs_every_particle_once():
    gap = fate_table([(1, 1, 3.0, 11.0, 0.0, 500.0, 8), (1, 3, 3.0, 11.0, 0.0, 500.0, 8)])
    with pytest.raises(ValueError, match="not 1"):
        release_fates(gap, classify_fates(gap, DOMAIN), [1], LAT, LON)
    missing = fate_table([(1, 1, 3.0, 11.0, 0.0, 500.0, 8)])
    with pytest.raises(ValueError, match="No fate records"):
        release_fates(missing, classify_fates(missing, DOMAIN), [1, 2], LAT, LON)


def small_day(ids, times, seed=0):
    """A tiny ``process_day``-shaped Dataset (only what ``assemble`` and the writer need)."""
    rng = np.random.default_rng(seed)
    n = len(ids)
    heights = (HEIGHT_BIN_EDGES_M[:-1] + HEIGHT_BIN_EDGES_M[1:]) / 2.0
    srr = rng.lognormal(-9, 2, size=(n, len(LAT), len(LON)))
    srr[:, 0, :] = 0.0
    data = {
        "srr": (("footprint", "latitude", "longitude"), srr),
        "particle_locations_top": (("footprint", "latitude", "longitude"), srr * 0.1),
        "mean_age_particles_top": (("footprint", "latitude", "longitude"), srr * 0 + 5.0),
        "particle_fraction_inside": ("footprint", np.full(n, 0.01)),
        "footprint_sounding_id": ("footprint", np.array(ids)),
        "time": ("sounding_id", pd.to_datetime(times)),
        "release_group": ("sounding_id", np.array(ids)),
        "air_temperature": ("sounding_id", np.full(n, 290.123456)),
        "column_weight": (("sounding_id", "release_level"), np.full((n, 2), 0.4)),
    }
    for face, along in (("n", "longitude"), ("e", "latitude"), ("s", "longitude"),
                        ("w", "latitude")):  # fmt: skip
        size = len(LON) if along == "longitude" else len(LAT)
        data[f"particle_locations_{face}"] = (
            ("footprint", "height", along),
            np.zeros((n, 20, size)),
        )
        data[f"mean_age_particles_{face}"] = (
            ("footprint", "height", along),
            np.zeros((n, 20, size)),
        )
    attrs = {"product": "GOSAT-1 v10.0", "site": "GOSAT-BRAZIL", "domain": "TEST",
             "l2_files": f"file_{times[0][:10]}", "max_level": 2, "particles_per_second": "30, 30",
             "model_version": "NAME III (version 8.7.1), ACRG fate-record patch (X20)",
             "met": "UMG_Mk6PT", "max_travel_hours": 720, "integration_window": "all_releases",
             "random_seed": "Fixed"}  # fmt: skip
    coords = {"sounding_id": np.array(ids), "footprint": np.arange(n), "latitude": LAT,
              "longitude": LON, "height": heights, "release_level": [1, 2]}  # fmt: skip
    return xr.Dataset(data, coords=coords, attrs=attrs)


def test_assemble_sorts_by_time_then_id_and_keeps_footprints_with_soundings():
    day2 = small_day(["b2", "a2"], ["2013-01-02 10:00", "2013-01-02 10:00"], seed=1)
    day1 = small_day(["z1"], ["2013-01-01 09:00"], seed=2)
    month = assemble([day2, day1])
    assert month["sounding_id"].values.tolist() == ["z1", "a2", "b2"]
    assert month["footprint_sounding_id"].values.tolist() == ["z1", "a2", "b2"]
    assert month["footprint_index"].values.tolist() == [0, 1, 2]
    np.testing.assert_array_equal(
        month["srr"].isel(footprint=1).values, day2["srr"].isel(footprint=1).values
    )
    with pytest.raises(ValueError, match="more than once"):
        assemble([day1, day1])


def test_written_file_rounds_relative_not_absolute(tmp_path):
    month = assemble([small_day(["a", "b", "c"], ["2013-01-05 10:00"] * 3)])
    path = tmp_path / "fp.nc"
    write_footprint_file(month, path)
    written = xr.open_dataset(path).load()
    srr = month["srr"].values
    nonzero = srr != 0
    relative = np.abs(written["srr"].values[nonzero] / srr[nonzero] - 1)
    # BitRound keeps BITROUND_BITS significant bits: error at most half a unit in the last bit.
    assert relative.max() <= 2.0 ** -(BITROUND_BITS + 1)
    assert (written["srr"].values[nonzero] != 0).all()  # small values are never zeroed
    assert (written["srr"].values[~nonzero] == 0).all()
    assert written["air_temperature"].values[0] == 290.123456  # small variables exact
    assert written["srr"].attrs["units"] == "(mol/mol)/(mol/m2/s)"
    assert written["srr"].attrs["loss_lifetime_hrs"] == -9.0
    assert written["height"].attrs["bounds"] == "height_bounds"
    assert written.attrs["species"] == "inert"
    assert written["footprint_index"].dtype == np.int32
    with pytest.raises(FileExistsError):
        write_footprint_file(month, path)
