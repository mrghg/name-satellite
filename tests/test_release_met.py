"""Tests for name_satellite.release_met (met at the release points)."""

import pandas as pd
import pytest

from name_satellite.release_met import (
    MAX_PER_GROUP,
    MAX_SOUNDINGS,
    release_met_blocks,
    release_met_soundings,
)


def plan(n_soundings, n_levels=3):
    rows = []
    for p in range(1, n_soundings + 1):
        for level in range(1, n_levels + 1):
            rows.append(
                {
                    "point_number": p,
                    "level_number": level,
                    "lon": -50.0 + p,
                    "lat": -10.0,
                    "release_end": pd.Timestamp("2013-01-01 15:00:00") + pd.Timedelta(seconds=p),
                    "z_pa": 100_000.0 - 10_000 * level - p,
                }
            )
    return pd.DataFrame(rows)


def test_one_row_per_sounding_at_level_1():
    soundings = release_met_soundings(plan(3))
    assert soundings["sounding"].tolist() == ["S0001", "S0002", "S0003"]
    assert soundings["z_pa"].tolist() == [89_999.0, 89_998.0, 89_997.0]
    assert soundings["group"].tolist() == [1, 1, 1]


def test_groups_and_limit():
    soundings = release_met_soundings(plan(MAX_PER_GROUP + 1, n_levels=1))
    assert soundings["group"].value_counts().to_dict() == {1: MAX_PER_GROUP, 2: 1}
    with pytest.raises(ValueError, match="at most 199"):
        release_met_soundings(plan(MAX_SOUNDINGS + 1, n_levels=1))


def test_blocks():
    text = release_met_blocks(release_met_soundings(plan(2)), 0.352, 0.234)
    assert "S0002, Lat-Long, -48.0, -10.0" in text
    assert "TGrid_S0002, 1, 00:01, 01/01/2013 15:00:02" in text
    assert "89999.0\n89998.0" in text  # pressures sorted, surface first
    assert text.count(", ReleaseMet_S0001") == 6
    assert "S0001_WD, Wind Direction (degrees), , , Lat-Long, HGrid_RM1, ZGrid_RM1," in text
    assert "S0001_PS, Pressure (Pa), , , , HGrid_RM1, ZGrid_Ground, TGrid_S0001," in text


# A release-met output file as NAME 8.7.1 writes it (cut down: 2 locations, 2 pressures).
TINY_OUTPUT = """\
 NAME III (version 8.7.1)
 Run name:                   test
 X grid size:                          2
 Number of preliminary cols:            4
 Number of field cols:                4

 Fields:
 , , , , , , , ,
 , , , , S0001_T, S0001_T, S0001_PS, S0001_BL,
 , , , , Temperature (K), Temperature (K), Pressure (Pa), Boundary Layer Depth,
 , , , , , , , ,
 , , , , K, K, Pa, m,
 , , , , , , , ,
 , , , , No averaging, No averaging, No averaging, No averaging,
 , , , , No averaging, No averaging, No averaging, No averaging,
 , , , , No averaging, No averaging, No averaging, No averaging,
 , , , , No averaging, No averaging, No averaging, No averaging,
 , , , , , , , ,
 , , , , , , , ,
 , , , , , , , ,
 , , , , Z = 100097.1 Pa, Z = 99886.84 Pa, Z = 0 m agl, Z = 0 m agl,
 , , , , , , , ,
 T Index, Location Index, T, Location, , , , ,
 1, 1, 01/01/2013 14:58 UTC, S0001, 299.032, 298.8596, 101360.8, 597.222,
 1, 2, 01/01/2013 14:58 UTC, S0002, 298.6956, 298.5284, 101000.0, 600.0,
"""


def test_read_release_met_file(tmp_path):
    from name_satellite.release_met import read_release_met_file

    path = tmp_path / "ReleaseMet_S0001_C1.txt"
    path.write_text(TINY_OUTPUT)
    columns, values = read_release_met_file(path)
    assert columns["request"].tolist() == ["S0001_T", "S0001_T", "S0001_PS", "S0001_BL"]
    assert columns["z"].tolist()[1] == "Z = 99886.84 Pa"
    assert values.loc["S0002", 1] == 298.5284


def test_release_met_per_source_pressure_rule():
    from name_satellite.release_met import release_met_per_source

    releases = plan(2)  # z_pa: level 1 ≈ 90,000 Pa, level 2 ≈ 80,000 Pa, level 3 ≈ 70,000 Pa
    met = pd.DataFrame(
        {
            "point_number": [1, 2],
            "temperature_k": [300.0, 290.0],
            "wind_speed_m_s": [5.0, 6.0],
            "wind_direction_deg": [90.0, 180.0],
            "boundary_layer_depth_m": [500.0, 600.0],
            "surface_pressure_pa": [95_000.0, 85_000.0],  # sounding 2: level 1 underground
        }
    )
    table = release_met_per_source(releases, met)
    assert table["pressure_pa"].tolist() == [89_999.0, 79_999.0, 69_999.0,
                                             85_000.0, 79_998.0, 69_998.0]  # fmt: skip
    assert table["temperature_c"].iloc[0] == pytest.approx(26.85)
    assert table["temperature_c"].iloc[1:3].isna().all()


def test_surface_pressure_only_blocks():
    from name_satellite.surface_pressure import QUANTITIES

    text = release_met_blocks(release_met_soundings(plan(2)), 0.352, 0.234, QUANTITIES, 20000.0)
    assert "ZArray" not in text  # no release-level grid needed
    assert "ZGrid_Top, m agl, 1, 0.0, 20000.0" in text  # domain top (D21)
    assert "S0001_PT, Pressure (Pa), , , , HGrid_RM1, ZGrid_Top," in text
    assert text.count(", ReleaseMet_S0001") == 2
    assert "S0002_PS, Pressure (Pa), , , , HGrid_RM1, ZGrid_Ground, TGrid_S0002," in text
