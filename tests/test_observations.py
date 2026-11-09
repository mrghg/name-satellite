"""Tests for name_satellite.observations (workflow step 1: L2 file -> sounding table)."""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from name_satellite.observations import (
    bad_pressure_mask,
    latlon_mask,
    legacy_bin,
    level_order_mask,
    quality_mask,
    read_l2,
    select,
    to_legacy_obs,
)
from name_satellite.releases import footprint_box_degrees

LEVELS = np.array([990.0, 940, 870, 800, 730, 650, 580, 510, 440, 370, 300, 230, 160, 86, 83, 80,
                   50, 10, 1, 0.1], dtype=np.float32)  # fmt: skip


def write_l2(path, lats, lons, times, quality=None, pressure=None):
    """A small UoL-style L2 file (n soundings × 20 levels)."""
    n = len(lats)
    pressure = np.tile(LEVELS, (n, 1)) if pressure is None else pressure
    ds = xr.Dataset(
        {
            "exposure_id": ("n", np.array([f"20130101{i:014d}" for i in range(n)], dtype="S22")),
            "time": ("n", pd.to_datetime(times).values),
            "latitude": ("n", np.array(lats, dtype=np.float32)),
            "longitude": ("n", np.array(lons, dtype=np.float32)),
            "xch4": ("n", np.arange(n, dtype=np.float32) + 1800),
            "xch4_uncertainty": ("n", np.full(n, 10.0, dtype=np.float32)),
            "xch4_quality_flag": ("n", np.zeros(n, np.int8) if quality is None else quality),
            "retr_flag": ("n", np.zeros(n, np.int8)),
            "pressure_levels": (("n", "m"), pressure),
            "pressure_weight": (("n", "m"), np.full((n, 20), 0.05, dtype=np.float32)),
            "xch4_averaging_kernel": (("n", "m"), np.ones((n, 20), dtype=np.float32)),
            "ch4_profile_apriori": (("n", "m"), np.full((n, 20), 1800.0, dtype=np.float32)),
        },
        attrs={"title": "NCEO GOSAT OCPR CH4", "bias_correction": "test offset"},
    )
    ds.to_netcdf(path)
    return path


def test_read_identifies_product_and_minute_times(tmp_path):
    path = write_l2(tmp_path / "UoL-GHG-L2-CH4-GOSAT-OCPR-20130101-fv10.0.nc", [-10, -11],
                    [-50, -51], ["2013-01-01 15:00", "2013-01-01 15:00"])  # fmt: skip
    table = read_l2(path)
    assert table.attrs["product"] == "UoL GOSAT OCPR"
    assert table.attrs["product_version"] == "10.0"
    assert table.attrs["time_resolution"].startswith("whole minutes")
    assert list(table["exposure_id"].values) == ["2013010100000000000000", "2013010100000000000001"]
    assert table.sizes == {"sounding": 2, "level": 20}


def test_read_refuses_unknown_files(tmp_path):
    path = write_l2(tmp_path / "something_else.nc", [0], [0], ["2013-01-01"])
    with pytest.raises(ValueError, match="Not a recognised"):
        read_l2(path)


def test_masks(tmp_path):
    pressure = np.tile(LEVELS, (4, 1))
    pressure[1, 0] = np.nan  # undefined level (as xarray reads a -9999.99 fill value)
    pressure[2, 0] = 930.0  # surface above level 1 (high ground)
    path = write_l2(tmp_path / "UoL-GHG-L2-CH4-GOSAT-OCPR-20130101-fv10.0.nc",
                    [-10, -11, -12, 7.251], [-50, -51, -52, -40], ["2013-01-01"] * 4,
                    quality=np.array([0, 0, 0, 1], np.int8), pressure=pressure)  # fmt: skip
    table = read_l2(path)
    assert quality_mask(table).tolist() == [True, True, True, False]
    assert bad_pressure_mask(table).tolist() == [True, False, True, True]
    assert level_order_mask(table).tolist() == [True, False, False, True]
    # Lower bound included, upper bound excluded (legacy latlon_filter).
    assert latlon_mask(table, (-11, 7.251), (-60, -30)).tolist() == [True, True, False, False]


def test_legacy_bin(tmp_path):
    # Soundings 0 and 1 share a 0.234 x 0.352 bin; 2 is in another bin; 3 is outside the bounds.
    path = write_l2(tmp_path / "UoL-GHG-L2-CH4-GOSAT-OCPR-20130101-fv10.0.nc",
                    [-10.00, -10.05, -12.0, 5.0], [-50.00, -49.90, -52.0, -40.0],
                    ["2013-01-01 15:00:00", "2013-01-01 15:00:30", "2013-01-01 14:00:00",
                     "2013-01-01 16:00:00"])  # fmt: skip
    table = read_l2(path)
    binned = legacy_bin(table, (-13.0, 0.0), (-55.0, -45.0), (0.234, 0.352))
    assert binned.sizes["sounding"] == 2
    first, second = binned.isel(sounding=0), binned.isel(sounding=1)  # sorted by time
    assert str(first["exposure_id"].values) == "2013010100000000000002"
    assert str(second["exposure_id"].values) == "2013010100000000000000,2013010100000000000001"
    assert second["time"].values == np.datetime64("2013-01-01T15:00:15")
    assert float(second["xch4"]) == pytest.approx(1800.5)
    # sqrt(mean(uncertainty)^2 + std(xch4)^2) with std over [1800, 1801] = 0.5
    assert float(second["xch4_uncertainty"]) == pytest.approx(np.sqrt(10.0**2 + 0.5**2))
    # Box = footprint + spread of the positions (float32, as stored in the files).
    lats = np.array([-10.0, -10.05], dtype=np.float32)
    lons = np.array([-50.0, -49.9], dtype=np.float32)
    dlon, dlat = footprint_box_degrees(np.mean(lats), 10.5)
    assert float(second["dlat"]) == pytest.approx(dlat + (lats.max() - lats.min()), rel=1e-12)
    assert float(second["dlon"]) == pytest.approx(dlon + (lons.max() - lons.min()), rel=1e-12)
    assert int(second["n_soundings"]) == 2


def test_to_legacy_obs(tmp_path):
    path = write_l2(tmp_path / "UoL-GHG-L2-CH4-GOSAT-OCPR-20130101-fv10.0.nc", [-10, -10.05],
                    [-50, -49.9], ["2013-01-01 15:00", "2013-01-01 15:01"])  # fmt: skip
    table = read_l2(path)
    single = to_legacy_obs(table)
    assert set(single.data_vars) >= {"lat", "lon", "pressure_weights", "exposure_id"}
    assert single["exposure_id"].dims == ("time", "id")
    binned = legacy_bin(table, (-13.0, 0.0), (-55.0, -45.0), (0.234, 0.352))
    obs = to_legacy_obs(select(binned, [True]))
    assert obs["exposure_id"].values.tolist() == [
        [b"2013010100000000000000", b"2013010100000000000001"]
    ]


def test_to_legacy_obs_refuses_duplicate_times(tmp_path):
    path = write_l2(tmp_path / "UoL-GHG-L2-CH4-GOSAT-OCPR-20130101-fv10.0.nc", [-10, -11],
                    [-50, -51], ["2013-01-01 15:00", "2013-01-01 15:00"])  # fmt: skip
    with pytest.raises(ValueError, match="Duplicate sounding times"):
        to_legacy_obs(read_l2(path))


# ---------------------------------------------------------------------------------------------
# Against the legacy preprocessing (acrg.satellite.gosat) on the same L2 files
# ---------------------------------------------------------------------------------------------

LEGACY_LAT_BOUNDS = (-35.753, 7.251)  # dev/scripts/run_legacy_gosat.py
LEGACY_LON_BOUNDS = (-75.984, -32.793)
LEGACY_BIN_SIZE = (0.234, 0.352)
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


def assert_same_soundings(ours, legacy, extra=()):
    """Same exposure IDs; every compared variable equal (dlat/dlon: legacy stored float32)."""
    our_ids = list(ours["exposure_id"].values)
    legacy_ids = list(legacy["exposure_id"].values)
    assert sorted(our_ids) == sorted(legacy_ids)
    ours = ours.isel(sounding=[our_ids.index(i) for i in legacy_ids])
    for name in COMPARED:
        np.testing.assert_array_equal(ours[name].values, legacy[name].values, err_msg=name)
    np.testing.assert_array_equal(ours["time"].values, legacy["time"].values)
    for name in extra:
        np.testing.assert_allclose(ours[name].values, legacy[name].values, rtol=1e-6, err_msg=name)


@pytest.mark.refdata
@pytest.mark.parametrize("day", ["20130101", "20130107"])
def test_matches_legacy_preprocessing(refdata, day):
    table = read_l2(refdata(f"cci_v70_l2_{day}"))
    table = select(table, quality_mask(table) & bad_pressure_mask(table))

    unbinned = select(table, latlon_mask(table, LEGACY_LAT_BOUNDS, LEGACY_LON_BOUNDS))
    legacy_unbinned = xr.open_dataset(refdata(f"cci_v70_legacy_unbinned_{day}"))
    assert_same_soundings(unbinned, legacy_unbinned)

    binned = legacy_bin(table, LEGACY_LAT_BOUNDS, LEGACY_LON_BOUNDS, LEGACY_BIN_SIZE)
    legacy_binned = xr.open_dataset(refdata(f"cci_v70_legacy_binned_{day}"))
    assert_same_soundings(binned, legacy_binned, extra=("dlat", "dlon"))
    assert binned.sizes["sounding"] == {"20130101": 25, "20130107": 64}[day]
