"""Tests for the future method (name_satellite.future): grids, input blocks, end positions."""

import numpy as np
import pandas as pd
import pytest

from name_satellite.future.grids import (
    FINE,
    MEDIUM,
    box_offset,
    box_shape,
    global_centres,
    name_grid_line,
    offset_from_centres,
    snapped_box,
)
from name_satellite.future.processing import N_END_HEIGHTS, end_position_fractions, in_block
from name_satellite.future.setup import Options, dispersion_options_block, hhmm


def test_snapped_boxes_lie_on_the_2p5_degree_lines_and_fit_their_blocks():
    rng = np.random.default_rng(1)
    for lon, lat in zip(rng.uniform(-150, 150, 200), rng.uniform(-80, 80, 200), strict=True):
        for level in (FINE, MEDIUM):
            box = snapped_box(lon, lat, level.half_width_deg)
            for edge in (box.west, box.east, box.south, box.north):
                assert edge % 2.5 == pytest.approx(0) or edge % 2.5 == pytest.approx(2.5)
            assert box.west <= lon - level.half_width_deg
            assert box.east >= lon + level.half_width_deg
            assert max(box_shape(box, level.cell_deg)) <= level.block_cells


def test_box_clipped_at_pole_and_dateline_refused():
    box = snapped_box(0.0, 88.0, 5.0)
    assert box.north == 90.0
    with pytest.raises(ValueError, match="dateline"):
        snapped_box(178.0, 0.0, 5.0)


def test_offset_round_trip_through_the_name_grid_line():
    box = snapped_box(-34.618, -10.390, 5.0)
    line = name_grid_line("HGridF1", box, FINE.cell_deg)
    _, _, nx, ny, _dx, _, x0, y0 = [v.strip() for v in line.split(",")]
    lat = [float(y0) + FINE.cell_deg * j for j in range(int(ny))]
    lon = [float(x0) + FINE.cell_deg * i for i in range(int(nx))]
    assert offset_from_centres(lat, lon, FINE.cell_deg) == box_offset(box, FINE.cell_deg)
    glat, glon = global_centres(FINE.cell_deg)
    j0, i0 = box_offset(box, FINE.cell_deg)
    assert glat[j0] == pytest.approx(lat[0])
    assert glon[i0] == pytest.approx(lon[0])


def test_in_block_pads_with_zeros():
    values = np.ones((3, 100, 125))
    out = in_block(values, FINE)
    assert out.shape == (3, 125, 125)
    assert out.sum() == values.sum()
    assert out[:, 100:].sum() == 0


def test_end_positions_count_every_particle_once():
    particles = pd.DataFrame(
        {
            "lon": [-60.0, -60.0, 179.99, 10.0],
            "lat": [-10.0, -10.0, 89.9, 0.0],
            "z": [500.0, 19999.0, 20500.0, 20001.0],
            "fate_code": [16, 16, 16, 8],
        }
    )
    hist = end_position_fractions(particles)
    assert hist.shape == (N_END_HEIGHTS, 72, 144)
    assert hist.sum() == pytest.approx(1.0)
    assert hist[-1].sum() == pytest.approx(0.25)  # the top exit
    assert hist[19].sum() == pytest.approx(0.5)  # 19,999 m and 20,500 m (still inside, clipped)


def test_dispersion_options_follow_name_rules():
    class S:
        max_particles = 100
        sync_minutes = 6

    text = dispersion_options_block(S, Options())
    values = [v.strip() for v in text.splitlines()[2].split(",")]
    assert values[1] == "7000000"
    assert values[2] == "01:00"
    assert values[4] == "01:00"
    assert values[10] == "New"
    long_range = dispersion_options_block(S, Options(near_source_hours=0.0))
    assert long_range.splitlines()[2].split(",")[1].strip() == "2"
    with pytest.raises(ValueError, match="skew"):
        dispersion_options_block(S, Options(near_source_hours=0.5, skew_hours=1.0))
    assert hhmm(1.5) == "01:30"
    assert hhmm(240) == "240:00"


def test_implausible_surface_pressure_flags_ratios_outside_the_limits():
    from name_satellite.surface_pressure import implausible_surface_pressure

    table = pd.DataFrame({"ratio": [0.865, 0.95, 1.0, 1.09, 1.15]})
    assert implausible_surface_pressure(table).tolist() == [True, False, False, False, True]


def test_input_refuses_more_field_requests_than_name_allows():
    from name_satellite.future.setup import MAX_FIELD_REQUESTS_PER_RUN, input_text

    plan = pd.DataFrame({"release_group": ["G"] * (MAX_FIELD_REQUESTS_PER_RUN // 3 + 1)})
    with pytest.raises(ValueError, match="field requests"):
        input_text(plan, None, None, Options(layout="sources"))


def test_species_layout_refuses_more_soundings_than_particles_can_carry():
    from name_satellite.future.setup import MAX_SPECIES_PER_RUN, input_text

    n = MAX_SPECIES_PER_RUN + 1
    weights = pd.DataFrame({"sounding_id": [f"s{k}" for k in range(n)], "release_group": "G",
                            "level_number": 1, "weight": 0.1})  # fmt: skip
    with pytest.raises(ValueError, match="species"):
        input_text(pd.DataFrame({"release_group": ["G"]}), None, None, Options(), weights)
    with pytest.raises(ValueError, match="weights"):
        input_text(pd.DataFrame({"release_group": ["G"]}), None, None, Options())


def test_split_plan_keeps_groups_whole_and_balances_releases():
    from name_satellite.future.setup import split_plan

    sizes = [16, 16, 32, 16, 16, 48, 16, 16]  # releases per group
    groups = [f"g{k}" for k, n in enumerate(sizes) for _ in range(n)]
    names = [f"SourceID_{n}" for n in range(1, len(groups) + 1)]
    plan = pd.DataFrame({"release_group": groups, "source_name": names})
    parts = [split_plan(plan, k, 4) for k in range(1, 5)]
    assert [
        n for part in parts for n in part["day_source_name"]
    ] == names  # every release once, in order
    for part in parts:
        assert list(part["source_name"]) == [f"SourceID_{n}" for n in range(1, len(part) + 1)]
        whole = plan[plan["release_group"].isin(part["release_group"])]
        assert len(whole) == len(part)  # no group divided
    assert [len(p) for p in parts] == [32, 48, 64, 32]
    with pytest.raises(ValueError, match="Part 5 of 4"):
        split_plan(plan, 5, 4)


def test_species_layout_weights_each_release_by_its_soundings():
    from types import SimpleNamespace

    from name_satellite.future.setup import (
        species_names,
        species_sources_and_requests_block,
    )

    # Two release groups of two levels: group A shared by soundings s1 and s2, group B by s3.
    plan = pd.DataFrame({
        "source_name": [f"SourceID_{k}" for k in range(1, 5)],
        "release_group": ["A", "A", "B", "B"], "level_number": [1, 2, 1, 2],
        "release_start": pd.Timestamp("2013-01-07 15:01"),
        "release_end": pd.Timestamp("2013-01-07 15:00"),
        "particles_per_second": 30, "lon": -50.0, "lat": -10.0, "z_pa": 90000.0, "dlon": 0.1,
        "dlat": 0.1, "dz_pa": 3000.0,
    })  # fmt: skip
    weights = pd.DataFrame({
        "sounding_id": ["s1", "s1", "s2", "s2", "s3", "s3"],
        "release_group": ["A", "A", "A", "A", "B", "B"],
        "level_number": [1, 2, 1, 2, 1, 2], "weight": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
    })  # fmt: skip
    boxes = pd.DataFrame({"release_group": ["A", "B"], "fine_grid": ["HGridF1", "HGridF2"],
                          "medium_grid": ["HGridM1", "HGridM2"]})  # fmt: skip
    species = species_names(weights)
    assert species.to_dict() == {"s1": "SND1", "s2": "SND2", "s3": "SND3"}
    text = species_sources_and_requests_block(
        plan, boxes, weights, species, SimpleNamespace(sync_minutes=6), {"window_hours": 242}
    )
    lines = text.splitlines()
    assert lines[1].count("Source Strength") == 2  # at most two soundings share a group
    assert ", SND1 2.000000e-01 g/s, SND2 4.000000e-01 g/s, Cuboid" in lines[3]  # A, level 2
    assert ", SND3 5.000000e-01 g/s, , Cuboid" in lines[4]  # B, level 1
    requests = [line for line in lines if line.startswith("Req_")]
    assert len(requests) == 9  # three grids per sounding
    assert "Req_HGridF1_SND2, Mixing Ratio, SND2, , HGridF1," in text
    assert "Req_HGridC_SND3, Mixing Ratio, SND3, , HGridC," in text


def test_main_job_uses_two_threads_and_memory_by_met_version(tmp_path):
    from name_satellite.future.setup import write_main_job
    from name_satellite.met import MET_VERSIONS

    by_name = {v.name: v for v in MET_VERSIONS}
    for name, memory in (("UMG_Mk6PT", "16G"), ("UMG_Mk7PT", "16G"), ("UMG_Mk10PT", "48G")):
        folder = tmp_path / name
        folder.mkdir()
        write_main_job(folder, "test", Options(), by_name[name])
        job = (folder / "job.sbatch").read_text()
        assert f"#SBATCH --mem={memory}\n" in job
        assert "#SBATCH --cpus-per-task=2\n" in job
