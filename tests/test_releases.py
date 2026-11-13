"""Tests for name_satellite.releases (workflow step 2: plan releases).

The golden tests check that we regenerate the vertical layers of every NAME source in the legacy
GOSAT-BRAZIL January 2013 runs from the observation files alone.
"""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from name_satellite.releases import (
    PARTICLE_RATES,
    build_release_plan,
    layer_boundaries_hpa,
    particles_per_second,
    release_layers_pa,
)
from name_satellite.validation.legacy_files import read_name_input_sources

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


def test_scaling_to_model_surface_keeps_each_layers_share_of_the_column():
    # D19: every level scaled by model / retrieval surface pressure.
    from name_satellite.releases import layer_boundaries_hpa, scale_to_model_surface

    levels = np.array([960.0, 900.0, 800.0, 600.0, 300.0, 100.0])
    scaled = scale_to_model_surface(levels, 970.0)
    assert scaled[0] == pytest.approx(970.0)
    original = layer_boundaries_hpa(levels)
    boundaries = layer_boundaries_hpa(scaled)
    assert boundaries[0] == pytest.approx(970.0)  # layer 0 starts at the model ground
    share = -np.diff(original) / original[0]
    assert -np.diff(boundaries) / boundaries[0] == pytest.approx(share)
    with pytest.raises(ValueError, match="> 0"):
        scale_to_model_surface(levels, 0.0)


def test_plan_scaled_to_model_surface():

    obs = xr.Dataset(
        {
            "time": ("t", pd.to_datetime(["2013-01-01 15:00:00", "2013-01-01 15:01:00"])),
            "lon": ("t", [-50.0, -51.0]),
            "lat": ("t", [-10.0, -11.0]),
            "pressure_levels": (("t", "lev"), np.tile(np.linspace(1000.0, 50.0, 20), (2, 1))),
        }
    )
    plan = build_release_plan(obs, model_surface_pressure_hpa=np.array([1010.0, 990.0]))
    level_1 = plan[plan["level_number"] == 1]
    assert (level_1["z_pa"] + level_1["dz_pa"] / 2).tolist() == pytest.approx([101_000, 99_000])
    with pytest.raises(ValueError, match="legacy"):
        build_release_plan(obs, surface_pressure_hpa=[1.0, 1.0],
                           model_surface_pressure_hpa=[1.0, 1.0])  # fmt: skip


def test_check_modelled_levels():
    # D21: the modelled layers must lie inside the domain and carry enough column weight.
    from name_satellite.releases import check_modelled_levels

    levels = np.array([[1000.0, 800.0, 600.0, 400.0, 200.0]])
    good = xr.Dataset(
        {
            "pressure_weights": (("s", "lev"), [[0.3, 0.3, 0.3, 0.05, 0.05]]),
            "xch4_averaging_kernel": (("s", "lev"), [[1.0] * 5]),
        },
        coords={"sounding_id": ("s", ["A"])},
    )
    check_modelled_levels(good, levels, [400.0], n_levels=3)  # layer 3 tops at 500 hPa: inside
    with pytest.raises(ValueError, match="above the domain top"):
        check_modelled_levels(good, levels, [550.0], n_levels=3)
    with pytest.raises(ValueError, match="column weight"):
        check_modelled_levels(good, levels, [100.0], n_levels=2)  # 0.6 < 0.85


def test_particle_rates_per_level():
    """D22: per-level rates replace the legacy level_split rule; one rate per modelled level."""
    obs = xr.Dataset(
        {
            "time": ("time", pd.to_datetime(["2013-01-01 15:00:00"])),
            "lon": ("time", [-50.0]),
            "lat": ("time", [-10.0]),
            "pressure_levels": (("time", "lev"), [np.linspace(1000.0, 50.0, 20)]),
        }
    )
    rates = PARTICLE_RATES["GOSAT-1 v10.0"]
    plan = build_release_plan(obs, n_levels=16, particle_rates=rates)
    assert plan["particles_per_second"].tolist() == [30] * 13 + [10, 3, 3]
    legacy = build_release_plan(obs, n_levels=16)
    assert legacy["particles_per_second"].tolist() == [30] * 12 + [3] * 4
    with pytest.raises(ValueError, match="16 positive particle rates"):
        build_release_plan(obs, n_levels=16, particle_rates=rates[:15])
