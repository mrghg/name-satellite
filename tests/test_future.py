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


def test_coarse_grid_copies_hold_whole_groups_and_at_most_the_limit():
    from name_satellite.future.setup import MAX_RELEASES_PER_COARSE_GRID, group_boxes

    groups = [f"G{g}" for g in range(95) for _ in range(16)]  # 95 groups of 16 releases
    plan = pd.DataFrame({"release_group": groups, "lon": -60.0, "lat": -10.0})
    boxes = group_boxes(plan)
    per_grid = boxes.groupby("coarse_grid").size() * 16
    assert (per_grid <= MAX_RELEASES_PER_COARSE_GRID).all()
    assert per_grid.sum() == len(plan)
    assert list(dict.fromkeys(boxes["coarse_grid"])) == ["HGridC1", "HGridC2", "HGridC3", "HGridC4"]
