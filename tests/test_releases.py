"""Tests for name_satellite.releases (workflow step 2: plan releases).

The golden tests check that we regenerate the vertical layers of every NAME source in the legacy
GOSAT-BRAZIL January 2013 runs from the observation files alone.
"""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from name_satellite.legacy import read_name_input_sources
from name_satellite.releases import (
    layer_boundaries_hpa,
    particles_per_second,
    release_layers_pa,
)

# ---------------------------------------------------------------------------------------------
# Small hand-checkable cases
# ---------------------------------------------------------------------------------------------


def test_boundaries_simple_profile():
    # Levels 1000, 800, 600 hPa: boundaries surface, halfway points, symmetric top.
    boundaries = layer_boundaries_hpa([1000.0, 800.0, 600.0])
    np.testing.assert_allclose(boundaries, [1000.0, 900.0, 700.0, 500.0])


def test_layers_simple_profile():
    centre, thickness = release_layers_pa([1000.0, 800.0, 600.0])
    # Layer 0: 1000-900 hPa; layer 1: 900-700 hPa; layer 2: 700-500 hPa (all in Pa).
    np.testing.assert_allclose(centre, [95000.0, 80000.0, 60000.0])
    np.testing.assert_allclose(thickness, [10000.0, 20000.0, 20000.0])


def test_layers_fill_column_without_gaps():
    levels = [1013.0, 950.0, 800.0, 500.0, 200.0, 80.0, 50.0, 10.0]
    centre, thickness = release_layers_pa(levels)
    # Total thickness = surface pressure minus the top boundary.
    top = layer_boundaries_hpa(levels)[-1]
    assert thickness.sum() == pytest.approx((1013.0 - top) * 100.0)
    # Adjacent layers touch: top of layer i == bottom of layer i+1.
    np.testing.assert_allclose(centre[:-1] - thickness[:-1] / 2, centre[1:] + thickness[1:] / 2)


def test_top_layer_stops_at_zero_pressure():
    # GOSAT's top levels are 1 and 0.1 hPa: the symmetric extension would go below 0 hPa.
    boundaries = layer_boundaries_hpa([1000.0, 10.0, 1.0, 0.1])
    assert boundaries[-1] == 0.0
    assert boundaries[-2] == pytest.approx(0.55)


def test_surface_pressure_override_changes_only_bottom_two_layers():
    levels = [990.0, 900.0, 800.0, 700.0]
    centre, thickness = release_layers_pa(levels)
    centre2, thickness2 = release_layers_pa(levels, surface_pressure_hpa=1010.0)
    np.testing.assert_array_equal(centre[2:], centre2[2:])
    np.testing.assert_array_equal(thickness[2:], thickness2[2:])
    assert not np.isclose(centre[0], centre2[0])
    assert not np.isclose(centre[1], centre2[1])
    assert thickness2[0] == pytest.approx((1010.0 - 955.0) * 100.0)  # surface to (1010+900)/2
    assert levels[0] == 990.0  # input not modified


def test_rejects_non_decreasing_levels():
    with pytest.raises(ValueError, match="decrease"):
        release_layers_pa([1000.0, 1000.0, 800.0])
    with pytest.raises(ValueError, match="decrease"):
        release_layers_pa([800.0, 900.0, 1000.0])


def test_particles_per_second_legacy_rule():
    assert [particles_per_second(level) for level in (1, 12, 13, 17)] == [30, 30, 3, 3]
    with pytest.raises(ValueError, match="1-based"):
        particles_per_second(0)


# ---------------------------------------------------------------------------------------------
# Golden tests against the legacy NAME input files
# ---------------------------------------------------------------------------------------------

N_LEVELS = 17  # legacy max_name_level
RTOL = 1e-6  # legacy files store about 9-12 significant figures


def check_day_against_legacy(input_file, obs_file):
    """Compare our layers with every legacy source of one day; return the inferred surface p."""
    sources = read_name_input_sources(input_file)
    obs = xr.open_dataset(obs_file)
    n_soundings = obs.sizes["time"]
    assert len(sources) == n_soundings * N_LEVELS

    # Legacy sources are ordered sounding by sounding, levels 1..17 within each.
    z = sources["z"].to_numpy().reshape(n_soundings, N_LEVELS)
    dz = sources["dz"].to_numpy().reshape(n_soundings, N_LEVELS)
    levels_hpa = obs["pressure_levels"].to_numpy().astype(float)

    # Horizontal position and timing.
    np.testing.assert_allclose(sources["x"], np.repeat(obs["lon"].to_numpy(), N_LEVELS), atol=1e-9)
    np.testing.assert_allclose(sources["y"], np.repeat(obs["lat"].to_numpy(), N_LEVELS), atol=1e-9)
    assert (sources["start_time"] - sources["stop_time"] == pd.Timedelta(seconds=60)).all()
    obs_times = np.repeat(obs["time"].to_numpy(), N_LEVELS)
    assert (abs(sources["stop_time"].to_numpy() - obs_times) < np.timedelta64(1, "s")).all()

    # Particle rates follow the legacy level rule.
    level_numbers = np.tile(np.arange(1, N_LEVELS + 1), n_soundings)
    expected_rates = [particles_per_second(level) for level in level_numbers]
    assert sources["particles_per_second"].tolist() == expected_rates

    # Vertical layers. Legacy used a NAME-derived surface pressure (X9), which we don't have, so:
    # (a) levels 3-17 don't depend on it and must match using the retrieval's own levels;
    for i in range(n_soundings):
        centre, thickness = release_layers_pa(levels_hpa[i])
        np.testing.assert_allclose(z[i, 2:], centre[2:N_LEVELS], rtol=RTOL)
        np.testing.assert_allclose(dz[i, 2:], thickness[2:N_LEVELS], rtol=RTOL)
    # (b) with the surface pressure implied by the legacy bottom layer, levels 1-2 match too.
    #     The inferred surface pressure carries the rounding of the legacy file (~1e-3 Pa), which
    #     is a large *relative* error for very thin bottom layers (~200 Pa), so use an absolute
    #     tolerance of 0.01 Pa (1e-4 hPa) here.
    surface_pa = z[:, 0] + dz[:, 0] / 2
    for i in range(n_soundings):
        centre, thickness = release_layers_pa(levels_hpa[i], surface_pa[i] / 100.0)
        np.testing.assert_allclose(z[i], centre[:N_LEVELS], rtol=0, atol=0.01)
        np.testing.assert_allclose(dz[i], thickness[:N_LEVELS], rtol=0, atol=0.01)
    return surface_pa / 100.0


@pytest.mark.refdata
def test_golden_mr2019_20130101(refdata):
    check_day_against_legacy(refdata("mr2019_input"), refdata("obs"))


@pytest.mark.refdata
@pytest.mark.slow
def test_golden_conc2018_all_january(refdata):
    run = refdata("conc2018_input").parents[1]
    n_checked = 0
    for day in range(1, 32):
        date = f"201301{day:02d}"
        input_file = run / "Input_files" / f"BackRun_SOUTHAMERICA_GOSAT-BRAZIL_column_{date}.txt"
        obs_file = run / "Observations" / f"gosat-fts_gosat_{date}_ch4-column.nc"
        n_checked += len(check_day_against_legacy(input_file, obs_file))
    assert n_checked == 678
