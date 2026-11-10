"""Tests for the release plan, met choice and NAME input writer (workflow steps 2-3)."""

import itertools

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from name_satellite.legacy import read_name_input_sources
from name_satellite.met import (
    MET_VERSIONS,
    met_version_for_period,
    required_met_files,
)
from name_satellite.name_input import (
    DispersionDomain,
    OutputGrid,
    RunSettings,
    name_input_text,
    run_timing,
)
from name_satellite.releases import build_release_plan, footprint_box_degrees

# The legacy SOUTHAMERICA set-up (from the 2019 legacy input file).
GRID = OutputGrid(nx=190, ny=357, dx=0.352, dy=0.234, x_min=-91.329, y_min=-60.981)
DOMAIN = DispersionDomain(x_min=-96.329, x_max=-19.449, y_min=-64.981, y_max=27.557)


def settings(**changes):
    values = {
        "run_name": "SOUTHAMERICA_GOSAT-BRAZIL_column_20130101",
        "output_folder": "/scratch/out",
        "grid": GRID,
        "domain": DOMAIN,
        "met_version": met_version_for_period("2012-12-01", "2013-01-01"),
        "met_folder": "/scratch/met",
    }
    values.update(changes)
    return RunSettings(**values)


# ---------------------------------------------------------------------------------------------
# Release plan and footprint box
# ---------------------------------------------------------------------------------------------


def test_footprint_box_at_equator_and_high_latitude():
    dlon, dlat = footprint_box_degrees(0.0)
    assert dlat == pytest.approx(np.degrees(10.5 / 6371.0), rel=1e-6)
    assert dlon == pytest.approx(dlat, rel=1e-9)
    dlon60, _ = footprint_box_degrees(60.0)
    assert dlon60 == pytest.approx(2 * dlat, rel=1e-3)  # 1/cos(60°) = 2


@pytest.mark.refdata
def test_release_plan_matches_legacy_sources(refdata):
    obs = xr.open_dataset(refdata("obs"))
    legacy = read_name_input_sources(refdata("mr2019_input"))
    plan = build_release_plan(obs)
    assert plan["source_name"].tolist() == legacy["name"].tolist()
    assert (plan["release_end"] == legacy["stop_time"]).all()
    assert (plan["release_start"] == legacy["start_time"]).all()
    assert plan["particles_per_second"].tolist() == legacy["particles_per_second"].tolist()
    # Un-binned soundings use the plain 10.5 km box (binned ones were widened by the bin spread,
    # which is not in the observation file).
    unbinned = np.isclose(legacy["dy"], 0.0944287686215)
    assert unbinned.sum() == 204
    np.testing.assert_allclose(plan["dlat"][unbinned], legacy["dy"][unbinned], rtol=1e-11)
    np.testing.assert_allclose(plan["dlon"][unbinned], legacy["dx"][unbinned], rtol=1e-11)


# ---------------------------------------------------------------------------------------------
# Met
# ---------------------------------------------------------------------------------------------


def test_met_version_choice():
    assert met_version_for_period("2012-12-01", "2013-01-01").name == "UMG_Mk6PT"
    assert met_version_for_period("2022-05-04", "2022-06-01").name == "UMG_Mk11PT"
    assert met_version_for_period("2026-02-01", "2026-03-01").name == "UMG_Mk12PT"
    with pytest.raises(ValueError, match="No single met version"):
        met_version_for_period("2013-04-01", "2013-05-15")  # Mk6 -> Mk7 changeover


def test_met_versions_are_contiguous():
    for earlier, later in itertools.pairwise(MET_VERSIONS):
        assert pd.Timestamp(earlier.last_day) + pd.Timedelta(days=1) == pd.Timestamp(
            later.first_day
        )


@pytest.mark.refdata
def test_required_met_files_for_january_2013():
    version = met_version_for_period("2012-12-01", "2013-01-01")
    files = required_met_files(version, "2012-12-01T14:58", "2013-01-01T16:46")
    # 3-hourly, 14 partitions, with one extra time step at each end.
    assert files[0] == "MO201212010900.UMG_Mk6_L59PT1.pp"
    assert len(files) == 14 * len(pd.date_range("2012-12-01T09:00", "2013-01-01T21:00", freq="3h"))


# ---------------------------------------------------------------------------------------------
# Input writer
# ---------------------------------------------------------------------------------------------


@pytest.fixture
def plan_20130101(refdata):
    return build_release_plan(xr.open_dataset(refdata("obs")))


@pytest.mark.refdata
def test_timing_legacy_2019_reproduces_legacy_input(plan_20130101):
    timing = run_timing(plan_20130101, settings(integration="legacy_2019"))
    # From the 2019 legacy input file: TGridEndOfRun t0, Av Time / Max Travel Time 721:00.
    assert timing["window_start"] == pd.Timestamp("2012-12-02 14:58:50")
    assert timing["window_hours"] == 721
    assert timing["max_travel_hours"] == 721


@pytest.mark.refdata
def test_timing_covers_every_release(plan_20130101):
    timing = run_timing(plan_20130101, settings())
    window_end = timing["window_start"] + pd.Timedelta(hours=timing["window_hours"])
    assert window_end >= plan_20130101["release_start"].max()
    # Every particle can travel its full 720 h inside the window.
    oldest = plan_20130101["release_end"].min() - pd.Timedelta(hours=timing["max_travel_hours"])
    assert oldest >= timing["window_start"]
    assert timing["window_hours"] == 722  # 720 h + the day's 1.8 h span, rounded up


@pytest.mark.refdata
def test_written_sources_round_trip(plan_20130101, tmp_path):
    path = tmp_path / "input.txt"
    path.write_text(name_input_text(plan_20130101, settings()))
    sources = read_name_input_sources(path)
    assert sources["name"].tolist() == plan_20130101["source_name"].tolist()
    for ours, theirs in [("lon", "x"), ("lat", "y"), ("z_pa", "z"), ("dlon", "dx"),
                         ("dlat", "dy"), ("dz_pa", "dz")]:  # fmt: skip
        # Written with full precision (repr), so the values come back exactly.
        np.testing.assert_array_equal(sources[theirs].to_numpy(), plan_20130101[ours].to_numpy())
    assert (sources["stop_time"] == plan_20130101["release_end"]).all()


@pytest.mark.refdata
def test_input_text_contents(plan_20130101):
    text = name_input_text(plan_20130101, settings(threads=4))
    assert text.count("\nSourceID_") == 425
    assert text.count("\nReq_H1Z1_") == 425
    assert "Mixing Ratio" in text
    assert "Particles, Yes, No, Yes" in text  # fate output on
    assert "Use OpenMP?, Threads, Output Group Threads\nYes, 4, 1" in text
    assert "MetDefnUMG_Mk6_L59PTpp.txt" in text
    assert "%" not in text.replace("Progress (%)", "")  # all placeholders filled
    # Legacy minimum boundary-layer depth of 40 m in every met module instance.
    assert text.count("Global_Mk6_PT1, 40.0,") == 1
    assert "TGridEndOfRun, 1, 01:00, 02/12/2012 14:58:50" in text


def test_field_requests_split_into_output_groups():
    from name_satellite.name_input import MAX_FIELD_COLUMNS_PER_GROUP, fields_output_group

    # 2,720 field columns in one group failed; 2,000 works (dev/notes/name-behaviour.md).
    assert MAX_FIELD_COLUMNS_PER_GROUP <= 2000
    assert fields_output_group(0) == "Fields_grid1"
    assert fields_output_group(MAX_FIELD_COLUMNS_PER_GROUP - 1) == "Fields_grid1"
    assert fields_output_group(MAX_FIELD_COLUMNS_PER_GROUP) == "Fields_grid1_2"


def test_restore_met_files(tmp_path):
    import gzip

    from name_satellite.met import MetVersion, restore_met_files

    archive = tmp_path / "archive"
    (archive / "Global" / "X").mkdir(parents=True)
    with gzip.open(archive / "Global" / "X" / "MO1.pp.gz", "wb") as f:
        f.write(b"zipped met")
    (archive / "Global" / "X" / "MO2.pp").write_bytes(b"plain met")
    version = MetVersion("X", "2000-01-01", None, "d.txt", "u.txt", "Global/X")

    met = tmp_path / "met"
    counts = restore_met_files(version, ["MO1.pp", "MO2.pp"], met, archive=archive)
    assert counts == {"restored": 1, "linked": 1, "already_present": 0}
    assert (met / "MO1.pp").read_bytes() == b"zipped met"
    assert (met / "MO2.pp").is_symlink()
    assert not list(met.glob("*.partial"))
    # Second call: nothing to do.
    again = restore_met_files(version, ["MO1.pp", "MO2.pp"], met, archive=archive)
    assert again["already_present"] == 2
    # A wrong-sized existing file is an error, not overwritten.
    (met / "MO1.pp").write_bytes(b"truncated")
    with pytest.raises(ValueError, match="size differs"):
        restore_met_files(version, ["MO1.pp"], met, archive=archive)
    with pytest.raises(FileNotFoundError):
        restore_met_files(version, ["MO3.pp"], met, archive=archive)
