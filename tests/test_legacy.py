"""Tests for name_satellite.legacy (readers for legacy pipeline files)."""

import numpy as np
import pytest

from name_satellite.legacy import read_legacy_met_folder, read_name_input_sources

N_LEVELS = 17


def check_met_matches_sources(met_folder, input_file):
    """Release-point met files line up one-to-one with the NAME sources (by point and level)."""
    met = read_legacy_met_folder(met_folder)
    sources = read_name_input_sources(input_file)
    assert len(met) == len(sources)
    # Source number = (point - 1) * 17 + level, i.e. sources are ordered point by point.
    source_number = (met["point_number"] - 1) * N_LEVELS + met["level_number"]
    assert source_number.tolist() == list(range(1, len(sources) + 1))
    # NAME prints 5 decimal places for positions and 7 significant figures for heights.
    np.testing.assert_allclose(met["lon"], sources["x"], atol=1e-5)
    np.testing.assert_allclose(met["lat"], sources["y"], atol=1e-5)
    np.testing.assert_allclose(met["z_pa"], sources["z"], rtol=1e-6)
    # Met is sampled at the release time (minute precision).
    assert (abs(met["time"] - sources["stop_time"]) < np.timedelta64(60, "s")).all()
    # Pressure at a release given in Pa is the requested pressure, EXCEPT where the requested
    # pressure is below the model's ground (higher than its surface pressure): then NAME samples
    # at the surface and reports a lower pressure. This happens for ~3% of legacy soundings
    # (docs/deviations.md X9). It is never higher.
    below_ground = met["pressure_pa"] < met["z_pa"] * (1 - 1e-6)
    np.testing.assert_allclose(
        met.loc[~below_ground, "pressure_pa"], met.loc[~below_ground, "z_pa"], rtol=1e-6
    )
    assert (met["pressure_pa"] <= met["z_pa"] * (1 + 1e-6)).all()
    assert met.loc[below_ground, "level_number"].isin([1, 2]).all()
    return int(below_ground.sum())


@pytest.mark.refdata
def test_mr2019_met(refdata):
    run = refdata("mr2019_input").parents[1]
    check_met_matches_sources(
        run / "Met" / "Met_SOUTHAMERICA_GOSAT-BRAZIL_column_20130101", refdata("mr2019_input")
    )


@pytest.mark.refdata
@pytest.mark.slow
def test_conc2018_met_all_january(refdata):
    run = refdata("conc2018_input").parents[1]
    n_below_ground = 0
    for day in range(1, 32):
        tag = f"SOUTHAMERICA_GOSAT-BRAZIL_column_201301{day:02d}"
        n_below_ground += check_met_matches_sources(
            run / "Met" / f"Met_{tag}", run / "Input_files" / f"BackRun_{tag}.txt"
        )
    assert n_below_ground == 22  # of 11,526 releases (recorded in docs/deviations.md X9)
