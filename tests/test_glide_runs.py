"""Tests for the GLIDE reference runs: spec numbers, release plan and NAME input options."""

import json
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from name_satellite.glide_runs import spec
from name_satellite.glide_runs.plan import point_release_plan
from name_satellite.glide_runs.prepare import prepare_run
from name_satellite.name_input import (
    ParticlePositionOutput,
    dispersion_options_block,
    met_declarations_block,
    name_input_text,
    particle_output_block,
    run_timing,
)

# ---------------------------------------------------------------------------------------------
# The spec's numbers (GLIDE docs/reference-runs.md)
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("domain", "n_cells", "first_centre", "last_centre"),
    [
        # §2 table: output cells, first and last cell centres (lon, lat)
        (spec.EUR, (279, 159), (-29.75, 32.25), (39.75, 71.75)),
        (spec.NAM, (259, 143), (-119.75, 22.25), (-55.25, 57.75)),
    ],
)
def test_output_grids_match_spec(domain, n_cells, first_centre, last_centre):
    g = domain.grid
    assert (g.nx, g.ny) == n_cells
    assert (g.dx, g.dy) == (0.25, 0.25)
    assert (g.x_min, g.y_min) == first_centre
    assert g.x_min + (g.nx - 1) * g.dx == pytest.approx(last_centre[0])
    assert g.y_min + (g.ny - 1) * g.dy == pytest.approx(last_centre[1])
    # The grid stops half a cell inside the particle domain (§2): its edges are 0.125 degrees
    # inside the box on every side.
    d = domain.dispersion
    assert g.x_min - g.dx / 2 == pytest.approx(d.x_min + 0.125)
    assert g.y_min - g.dy / 2 == pytest.approx(d.y_min + 0.125)
    assert g.x_min + (g.nx - 0.5) * g.dx == pytest.approx(d.x_max - 0.125)
    assert g.y_min + (g.ny - 0.5) * g.dy == pytest.approx(d.y_max - 0.125)
    assert d.z_max_m_agl >= 20000.0  # no particle removed below 20 km above the ground (§2)


def test_ERA5_grid_points_are_cell_centres():
    # Cell centres sit on ERA5's 0.25-degree grid (§2).
    for domain in spec.DOMAINS.values():
        g = domain.grid
        assert (g.x_min / 0.25) % 1 == 0
        assert (g.y_min / 0.25) % 1 == 0


def test_blocks_and_release_count():
    for block in spec.BLOCKS.values():
        times = block.release_starts()
        assert len(times) == 336  # §4
        assert times[0] >= block.cover_start + pd.Timedelta(hours=spec.MAX_AGE_HOURS)
        assert times[-1] + pd.Timedelta(hours=1) <= block.cover_end
    receptors_per_domain = {
        name: sum(r.domain == name for r in spec.RECEPTORS.values()) for name in spec.DOMAINS
    }
    total = sum(
        receptors_per_domain[block.domain] * len(block.release_starts())
        for block in spec.BLOCKS.values()
    )
    assert total == 4704  # §4


def test_particle_rates():
    # 21,600 per point release and 100,800 per column release (§4).
    assert spec.POINT_PARTICLES_PER_SECOND * 3600 == 21600
    assert spec.COLUMN_PARTICLES_PER_SECOND * 3600 == 100800


def test_variants_change_only_their_settings():
    settings = spec.reference_settings("r", "/out", spec.EUR, "/met", spec.MAX_PARTICLES_POINT)
    assert settings.deep_convection == "New"
    assert (settings.min_bl_depth_m, settings.max_bl_depth_m) == (100.0, 4500.0)
    assert settings.sync_minutes == 5
    assert settings.random_seed == "Fixed"
    group = spec.variant_settings(settings, "group")
    assert group.deep_convection == "No"
    assert group.sync_minutes == 3
    assert (group.min_bl_depth_m, group.max_bl_depth_m) == (40.0, 4000.0)
    assert replace(group, deep_convection="New", sync_minutes=5, min_bl_depth_m=100.0,
                   max_bl_depth_m=4500.0) == settings  # fmt: skip
    assert spec.variant_settings(settings, "seed3").random_seed == "Random"
    assert spec.variant_settings(settings, "nomeso").mesoscale_motions is False
    with pytest.raises(ValueError, match="Unknown run variant"):
        spec.variant_settings(settings, "kernel")  # a FLEXPART-only run


# ---------------------------------------------------------------------------------------------
# Release plan
# ---------------------------------------------------------------------------------------------


def test_point_release_plan():
    times = pd.date_range("2019-01-08 00:00", "2019-01-08 23:00", freq="1h")
    plan = point_release_plan(spec.RECEPTORS["JFJ"], times)
    assert len(plan) == 24
    assert plan["source_name"].iloc[0] == "SourceID_1"
    # Backwards run: NAME's Start Time is the end of the hour named by release_time.
    assert (plan["release_start"] - plan["release_time"] == pd.Timedelta(hours=1)).all()
    assert (plan["release_end"] == plan["release_time"]).all()
    # A point 1,000 m above the model's ground, not 3,580 m above sea level (§3).
    assert (plan["z_coord"] == "m agl").all()
    assert (plan["z"] == 1000.0).all()
    assert (plan[["dz", "dlon", "dlat"]].to_numpy() == 0).all()


def test_point_release_plan_rejects_bad_input():
    mhd = spec.RECEPTORS["MHD"]
    with pytest.raises(ValueError, match="column receptor"):
        point_release_plan(spec.RECEPTORS["SGP"], ["2019-07-08 19:00"])  # column (§5)
    with pytest.raises(ValueError, match="on the hour"):
        point_release_plan(mhd, ["2019-01-08 00:30"])  # not on the hour
    with pytest.raises(ValueError, match="increasing"):
        point_release_plan(mhd, ["2019-01-08 01:00", "2019-01-08 00:00"])  # not increasing


def test_block_timing_covers_240_hours_from_every_release():
    block = spec.BLOCKS["EUR-Jan"]
    plan = point_release_plan(spec.RECEPTORS["MHD"], block.release_starts())
    settings = spec.reference_settings("r", "/out", spec.EUR, "/met", spec.MAX_PARTICLES_POINT)
    timing = run_timing(plan, settings)
    # Integration window: 240 h before the first release to the end of the last.
    assert timing["window_start"] == pd.Timestamp("2018-12-22 00:00") == block.cover_start
    assert timing["window_hours"] == 240 + 336
    assert timing["max_travel_hours"] == 240


# ---------------------------------------------------------------------------------------------
# NAME input options added for these runs (the satellite defaults must not change)
# ---------------------------------------------------------------------------------------------


def _settings(**changes):
    settings = spec.reference_settings("r", "/out", spec.EUR, "/met", spec.MAX_PARTICLES_POINT)
    return replace(settings, **changes)


def test_dispersion_options():
    line = dispersion_options_block(_settings()).splitlines()[-1]
    assert "00:05, 1, , Dispersion Domain, New, No, No, No, No, Yes, Yes, No" in line
    line = dispersion_options_block(_settings(mesoscale_motions=False)).splitlines()[-1]
    assert line.endswith("Yes, No, No")
    with pytest.raises(ValueError, match="Deep convection"):
        dispersion_options_block(_settings(deep_convection="Yes"))  # NAME wants No/Old/New


def test_satellite_dispersion_options_unchanged():
    from name_satellite.domains import SOUTHAMERICA
    from name_satellite.name_input import RunSettings

    settings = RunSettings(
        run_name="r", output_folder="/o", grid=SOUTHAMERICA.grid,
        domain=SOUTHAMERICA.dispersion, met_version=spec.MET_VERSION, met_folder="/m",
    )  # fmt: skip
    line = dispersion_options_block(settings).splitlines()[-1]
    assert line.endswith("00:06, 1, , Dispersion Domain, No, No, No, No, No, Yes, Yes, No")
    block = particle_output_block(settings)
    assert "T-Grid" not in block
    assert "Positions" not in block


def test_boundary_layer_limits_in_met_declarations():
    text = met_declarations_block(_settings())
    rows = [line for line in text.splitlines() if line.startswith("Global_Mk10_PT")]
    rows = [r for r in rows if "Flow" not in r]
    assert rows
    assert all([c.strip() for c in r.split(",")[1:3]] == ["100.0", "4500.0"] for r in rows)


def test_point_sources_and_positions_in_input():
    times = pd.date_range("2019-01-08 00:00", periods=3, freq="1h")
    plan = point_release_plan(spec.RECEPTORS["MHD"], times)
    positions = ParticlePositionOutput(pd.Timestamp("2019-01-08 03:00"), 4, 60)
    text = name_input_text(plan, _settings(max_age_hours=1, particle_positions=positions))
    source = next(line for line in text.splitlines() if line.startswith("SourceID_2,"))
    assert source.startswith(
        "SourceID_2, 08/01/2019 02:00:00, 08/01/2019 01:00:00, 6/sec, INERT-TRACER 1.0 g/s, "
        "Cuboid, Lat-Long, m agl, , , No, Yes, -9.9046, 53.3267, 10.0, 0.0, 0.0, 0.0,"
    )
    assert "TGridPositions, 4, 01:00, 08/01/2019 03:00:00" in text
    assert "Positions, Yes, No, No, No, No, No, No, No, Lat-Long, m agl, TGridPositions" in text


def test_prepare_run_writes_folder(tmp_path):
    folder = prepare_run(
        tmp_path / "run", "CBW-200", "2019-01-08T00", "2019-01-08T01", max_age_hours=1, qos="debug"
    )
    for name in ("input.txt", "job.sbatch", "plan.csv", "restore_met.sbatch", "run.json"):
        assert (folder / name).exists()
    text = (folder / "input.txt").read_text()
    assert "R1_TOPO, Topography" in text
    assert "R1_LAND, Land Fraction" in text
    assert "NAME-8.7.1-fateid" in (folder / "job.sbatch").read_text()
    record = json.loads((folder / "run.json").read_text())
    assert record["receptor"]["height_m_agl"] == 200.0
    assert record["n_releases"] == 2
    with pytest.raises(FileExistsError):
        prepare_run(folder, "CBW-200", "2019-01-08T00", "2019-01-08T01")
    plan = pd.read_csv(folder / "plan.csv")
    np.testing.assert_array_equal(plan["z"], [200.0, 200.0])


# ---------------------------------------------------------------------------------------------
# The Lamont column (§5)
# ---------------------------------------------------------------------------------------------


def _surface_pressure_table():
    return pd.DataFrame(
        {
            "release_time": pd.date_range("2019-07-08 18:00", periods=2, freq="1h"),
            "surface_pressure_pa": [97700.0, 97643.39],
            "domain_top_pressure_pa": [5500.0, 5500.0],
        }
    )


def test_column_release_plan():
    from name_satellite.glide_runs.column import column_release_plan

    plan = column_release_plan(spec.RECEPTORS["SGP"], _surface_pressure_table())
    # One source per release, from the model surface pressure to 100 hPa, even in pressure.
    assert (plan["z_coord"] == "Pa").all()
    np.testing.assert_allclose(plan["z"] - plan["dz"] / 2, [10000.0, 10000.0])
    np.testing.assert_allclose(plan["z"] + plan["dz"] / 2, [97700.0, 97643.39])
    assert (plan["particles_per_second"] == 28).all()
    with pytest.raises(ValueError, match="not a column"):
        column_release_plan(spec.RECEPTORS["LEF"], _surface_pressure_table())


def test_column_run_needs_surface_pressure(tmp_path):
    with pytest.raises(ValueError, match="is a column"):
        prepare_run(tmp_path / "a", "SGP", "2019-07-08T18", "2019-07-08T19")
    csv = tmp_path / "ps.csv"
    _surface_pressure_table().to_csv(csv, index=False)
    with pytest.raises(ValueError, match="lacks some"):
        prepare_run(
            tmp_path / "b", "SGP", "2019-07-08T18", "2019-07-08T20", surface_pressure_csv=csv
        )
    folder = prepare_run(
        tmp_path / "c", "SGP", "2019-07-08T18", "2019-07-08T19", surface_pressure_csv=csv
    )
    text = (folder / "input.txt").read_text()
    assert "28/sec, INERT-TRACER 1.0 g/s, Cuboid, Lat-Long, Pa," in text
    assert "6/sec" not in text
    assert "column_method" in json.loads((folder / "run.json").read_text())
