"""Tests for workflow step 6 (footprints, boundary histograms, column combination) and the legacy
reproduction of the processed footprint files."""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from name_satellite.boundary_histograms import (
    LEGACY_HEIGHT_CENTRES_M,
    _legacy_copy_source,
    legacy_boundary_histograms,
)
from name_satellite.column import combine_levels, nearest_release_levels
from name_satellite.footprints import (
    EARTH_RADIUS_M,
    footprint_from_mixing_ratio,
    grid_cell_area_m2,
)
from name_satellite.legacy_pipeline import process_day_legacy, write_legacy_fp_file

# ---------------------------------------------------------------------------------------------
# Footprints
# ---------------------------------------------------------------------------------------------


def test_cell_areas_cover_the_globe():
    lat = np.arange(-89.5, 90, 1.0)
    lon = np.arange(0.5, 360, 1.0)
    area = grid_cell_area_m2(lat, lon)
    assert area.sum() == pytest.approx(4 * np.pi * EARTH_RADIUS_M**2, rel=1e-12)
    # Cells shrink towards the poles.
    assert area[0, 0] < area[90, 0]


def test_footprint_conversion():
    # 1 ppm s over a 1e6 m2 cell from a 60 s release: 1e-6 * 1e6 / 60 (mol/mol)/(mol/m2/s).
    field = np.array([[1.0, 2.0]])
    area = np.array([[1e6, 1e6]])
    np.testing.assert_allclose(footprint_from_mixing_ratio(field, area, 60.0), [[1 / 60, 2 / 60]])


# ---------------------------------------------------------------------------------------------
# Boundary histograms (legacy)
# ---------------------------------------------------------------------------------------------

# A 4 x 3 output grid: lat 0, 1, 2, 3 and lon 10, 11, 12 (cell centres, 1° cells).
LAT = np.array([0.0, 1.0, 2.0, 3.0])
LON = np.array([10.0, 11.0, 12.0])


def particles(rows):
    return pd.DataFrame(rows, columns=["lon", "lat", "height_m_agl", "source_number"])


def test_histograms_simple_edges():
    p = particles(
        [
            [11.0, 3.7, 500.0, 1],  # north edge, middle lon cell, 0-1 km
            [12.9, 1.0, 2500.0, 1],  # east edge, lat cell 1, 2-3 km
            [11.0, 1.5, 3000.0, 1],  # inside the domain: not counted
        ]
    )
    h = legacy_boundary_histograms(p, LAT, LON, n_sources=1)
    assert h["north"][0, 1, 0] == 0.5
    assert h["east"][0, 1, 2] == 0.5
    assert h["south"].sum() == 0
    assert h["west"].sum() == 0


def test_histograms_legacy_corner_and_overshoot_quirks():
    p = particles(
        [
            [13.2, 3.8, 500.0, 1],  # beyond the NE corner: counted on BOTH north and east (X14)
            [9.0, 3.8, 500.0, 1],  # north edge but west of the first lon cell: north drops it
            [11.0, 0.0, 500.0, 1],  # ordinary south-edge particle
        ]
    )
    h19 = legacy_boundary_histograms(p, LAT, LON, n_sources=1, variant="2019")
    # 2019: bins end at the grid edge (lon 12.5, lat 3.5). The NE particle (13.2, 3.8) is beyond the
    # last bin of both the north edge (lon) and the east edge (lat), so it is dropped from both.
    # The particle at lon 9.0 is on the north edge (dropped: west of the first lon bin) and the
    # west edge (dropped: lat 3.8 beyond the last lat bin). Only the south particle is counted.
    assert h19["south"][0, 1, 0] == 1.0
    assert sum(h19[e].sum() for e in ("north", "east", "south", "west")) == pytest.approx(1.0)

    h25 = legacy_boundary_histograms(p, LAT, LON, n_sources=1, variant="2025")
    # 2025: bins stretch at the FAR end, so the NE particle now counts on north AND east, and the
    # west-edge particle (far end of the lat bins) counts on west; the north edge still drops the
    # particle west of the first lon cell (near end).
    assert h25["north"][0, 2, 0] == pytest.approx(1 / 4)
    assert h25["east"][0, 3, 0] == pytest.approx(1 / 4)
    assert h25["west"][0, 3, 0] == pytest.approx(1 / 4)
    assert h25["south"][0, 1, 0] == pytest.approx(1 / 4)


def test_histograms_copy_rule_when_no_particle_reaches_an_edge():
    n = 17 * 2  # two soundings
    rows = [[11.0, 0.0, 500.0, s] for s in range(1, n + 1) if s not in (5, 18)]
    h = legacy_boundary_histograms(particles(rows), LAT, LON, n_sources=n)
    assert h["copied_from"][4] == 4  # level 5 copies level 4 of the same sounding
    assert h["copied_from"][17] == 1  # sounding 2, level 1 copies sounding 1, level 1
    assert h["copied_from"].sum() == 5
    np.testing.assert_array_equal(h["south"][4], h["south"][3])


def test_legacy_copy_source_rule():
    assert _legacy_copy_source(2) == 1
    assert _legacy_copy_source(18) == 1  # level 1 of sounding 2 -> level 1 of sounding 1
    assert _legacy_copy_source(35) == 18  # level 1 of sounding 3 -> level 1 of sounding 2


def test_legacy_height_bins():
    assert LEGACY_HEIGHT_CENTRES_M[0] == 500
    assert LEGACY_HEIGHT_CENTRES_M[-1] == 19500
    assert len(LEGACY_HEIGHT_CENTRES_M) == 20


# ---------------------------------------------------------------------------------------------
# Column combination
# ---------------------------------------------------------------------------------------------


def test_nearest_release_levels():
    release_pa = np.array([99000.0, 90000.0, 80000.0])
    retrieval_hpa = np.array([1000.0, 880.0, 860.0, 700.0])
    # 860 hPa is nearer 900 hPa (40 away) than 800 hPa (60 away), so the 880 and 860 hPa retrieval
    # levels both use NAME release 1: the legacy nearest-pressure rule can use a release twice.
    np.testing.assert_array_equal(
        nearest_release_levels(release_pa, retrieval_hpa, 4), [0, 1, 1, 2]
    )
    np.testing.assert_array_equal(nearest_release_levels(release_pa, retrieval_hpa, 2), [0, 1])


def test_combine_levels():
    per_level = np.array([[1.0, 0.0], [0.0, 1.0], [2.0, 2.0]])
    result = combine_levels(per_level, [0, 2, 2], [0.5, 0.25, 0.25], [1.0, 1.0, 2.0])
    np.testing.assert_allclose(result, [0.5 + 0.5 + 1.0, 0.0 + 0.5 + 1.0])


# ---------------------------------------------------------------------------------------------
# Golden test: rebuild the 2019 production footprints for 1 Jan 2013 bit for bit
# ---------------------------------------------------------------------------------------------

COMPARED_VARIABLES = [
    "fp",
    "temperature",
    "pressure",
    "wind_speed",
    "wind_direction",
    "PBLH",
    "release_lon",
    "release_lat",
    "particle_locations_n",
    "particle_locations_e",
    "particle_locations_s",
    "particle_locations_w",
    "lat",
    "lon",
    "height",
]


@pytest.mark.refdata
def test_golden_reproduces_production_20130101(refdata, tmp_path):
    run = refdata("mr2019_input").parents[1]
    result = process_day_legacy(run, "SOUTHAMERICA_GOSAT-BRAZIL_column", "20130101")
    out = tmp_path / "fp_20130101.nc"
    write_legacy_fp_file(result, out)

    new = xr.open_dataset(out)
    production = xr.open_dataset(refdata("production_201301"))
    production = production.sel(time=new["time"].values)
    assert new.sizes["time"] == 25
    for name in COMPARED_VARIABLES:
        np.testing.assert_array_equal(new[name].values, production[name].values, err_msg=name)
    assert new.attrs["max_level"] == production.attrs["max_level"]


# Days of January 2013 on which ONE footprint cell differs from production by one rounding step of
# the legacy file writer (least_significant_digit=5 rounds to multiples of 2**-17; X7). Each such
# value lies within ~1e-8 of a rounding boundary, so the direction of rounding depends on the last
# bits of the floating-point arithmetic (and possibly on the 2019 numpy/netCDF4 versions). Neither
# rounding the output footprint nor the input fields to float32 (as legacy held the fields)
# reproduces all 31 days; each fixes some days and breaks others.
ROUNDING_STEP_DAYS = {"20130103", "20130110", "20130111", "20130113", "20130125"}


@pytest.mark.refdata
@pytest.mark.slow
def test_golden_month_matches_production(refdata, tmp_path):
    """Legacy processing of all 31 days of the 2019 MR run reproduces the production file:
    every variable exactly, except one footprint cell (one rounding step) on ROUNDING_STEP_DAYS."""
    from name_satellite.legacy import read_legacy_met_folder
    from name_satellite.legacy_pipeline import process_legacy_inputs
    from name_satellite.name_output import read_acrg_particle_file, read_fields

    run = refdata("mr2019_month_fields_20130131").parents[1]
    tag = "SOUTHAMERICA_GOSAT-BRAZIL_column"
    production = xr.open_dataset(refdata("production_201301")).load()
    step = 2.0**-17
    days_with_differences = set()
    for d in range(1, 32):
        day = f"201301{d:02d}"
        result = process_legacy_inputs(
            read_fields(run / "MixR_files" / f"Fields_{tag}_{day}.txt.gz"),
            read_acrg_particle_file(
                run / "Particle_files" / f"particle_location_{tag}_{day}.txt.gz"
            ),
            read_legacy_met_folder(run / "Met" / f"Met_{tag}_{day}"),
            xr.open_dataset(run / "Observations" / f"gosat-fts_gosat_{day}_ch4-column.nc").load(),
            release_seconds=60.0,  # legacy time_step.txt; not in the month archive
        )
        out = tmp_path / f"fp_{day}.nc"
        write_legacy_fp_file(result, out)
        new = xr.open_dataset(out).load()
        prod = production.sel(time=new["time"].values)
        for name in COMPARED_VARIABLES:
            if name == "fp":
                continue
            np.testing.assert_array_equal(
                new[name].values, prod[name].values, err_msg=f"{day} {name}"
            )
        difference = np.abs(new["fp"].values.astype(float) - prod["fp"].values.astype(float))
        n_cells = int((difference > 0).sum())
        if n_cells:
            days_with_differences.add(day)
            assert n_cells == 1, day
            assert difference.max() == pytest.approx(step, rel=1e-6), day
    assert days_with_differences == ROUNDING_STEP_DAYS
