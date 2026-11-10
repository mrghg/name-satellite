"""Tests for preparing NAME runs from L2 data: l2_files, domains, lotus jobs and prepare."""

import pandas as pd
import pytest
import xarray as xr
from test_observations import write_l2  # tests/ is on sys.path (no package)

from name_satellite import prepare
from name_satellite.domains import GOSAT_BRAZIL, SITES
from name_satellite.l2_files import L2FileMissing, l2_url
from name_satellite.lotus import write_job_script, write_met_restore_job


def test_l2_urls_and_coverage():
    assert l2_url("GOSAT-1 v10.0", "2013-01-07") == (
        "https://dap.ceda.ac.uk/neodc/gosat/data/ch4/UoL_gosat1_v10/CH4_GOSAT_OCPR/v10.0/2013/"
        "UoL-GHG-L2-CH4-GOSAT-OCPR-20130107-fv10.0.nc"
    )
    assert l2_url("GOSAT-2 v2.0", "2020-01-01").endswith(
        "UoL_gosat2_v2/CH4_GOSAT2_OCPR/v2.0/2020/UoL-GHG-L2-CH4-GOSAT2-OCPR-20200101-fv2.0.nc"
    )
    with pytest.raises(ValueError, match="outside GOSAT-2"):
        l2_url("GOSAT-2 v2.0", "2013-01-01")


def test_site_definitions():
    assert SITES["GOSAT-BRAZIL"] is GOSAT_BRAZIL
    grid = GOSAT_BRAZIL.domain.grid
    assert (grid.nx, grid.ny, grid.x_min, grid.y_min) == (190, 357, -91.329, -60.981)


def test_job_script_fails_on_name_errors(tmp_path):
    text = write_job_script(tmp_path / "job", tmp_path, tmp_path / "input.txt", threads=1)
    text = text.read_text()
    assert f"error_file={tmp_path / 'inputError.txt'}" in text
    assert text.rstrip().endswith("exit $status")
    restore = write_met_restore_job(tmp_path / "r", tmp_path, "UMG_Mk6PT", "a", "b", "/m")
    assert "python -m name_satellite.met_restore" in restore.read_text()


def test_prepare_day(tmp_path, monkeypatch):
    # Three soundings: two in the GOSAT-BRAZIL box (one sharing a time), one outside.
    l2 = write_l2(
        tmp_path / "UoL-GHG-L2-CH4-GOSAT-OCPR-20130101-fv10.0.nc",
        [-10.0, -11.0, 20.0],
        [-50.0, -51.0, -50.0],
        ["2013-01-01 15:00", "2013-01-01 15:00", "2013-01-01 15:01"],
    )
    monkeypatch.setattr(prepare, "download_l2", lambda product, day, cache: l2)
    row = prepare.prepare_day("GOSAT-1 v10.0", GOSAT_BRAZIL, "2013-01-01", tmp_path / "month")
    assert row["status"] == "prepared"
    assert (row["l2_soundings"], row["soundings"], row["sources"]) == (3, 2, 34)

    run = tmp_path / "month" / "20130101"
    obs = xr.open_dataset(run / "obs_GOSAT-BRAZIL_20130101.nc")
    assert obs.sizes["sounding_id"] == 2
    plan = pd.read_csv(run / "release_plan.csv")
    assert sorted(plan["sounding_id"].astype(str).unique()) == sorted(obs["sounding_id"].values)
    text = (run / "input.txt").read_text()
    assert text.count("\nSourceID_") == 34
    assert "ReleaseMet_S0002" in text  # release-point met requested in the run (X21)


def test_prepare_day_without_l2_file(tmp_path, monkeypatch):
    def missing(product, day, cache):
        raise L2FileMissing("none")

    monkeypatch.setattr(prepare, "download_l2", missing)
    row = prepare.prepare_day("GOSAT-1 v10.0", GOSAT_BRAZIL, "2013-01-06", tmp_path)
    assert row["status"] == "no L2 file"
