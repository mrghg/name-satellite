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


def surface_pressure_output(name, names, pressures_pa, top_pa=5000.0):
    """A met-only run's output file as NAME 8.7.1 writes it: surface pressure and the pressure at
    the domain top (20,000 m agl) at each location."""
    rows = [f" 1, {i}, 01/01/2013 15:00 UTC, {n}, {p}, {top_pa}," for i, (n, p) in
            enumerate(zip(names, pressures_pa, strict=True), start=1)]  # fmt: skip
    return "\n".join(
        [" NAME III (version 8.7.1)", " Fields:", " , , , , , ,", f" , , , , {name}_PS, {name}_PT,",
         " , , , , Pressure (Pa), Pressure (Pa),", " , , , , , ,", " , , , , Pa, Pa,",
         " , , , , Z = 0 m agl, Z = 20000 m agl,", " T Index, Location Index, T, Location, , ,",
         *rows, ""]
    )  # fmt: skip


# Pressure weights with 0.96 of the column in the 16 modelled levels (D21).
WEIGHTS = [0.06] * 16 + [0.01] * 4


def test_prepare_day_then_main_input(tmp_path, monkeypatch):
    # Three soundings: two in the GOSAT-BRAZIL box (sharing a time), one outside.
    l2 = write_l2(
        tmp_path / "UoL-GHG-L2-CH4-GOSAT-OCPR-20130101-fv10.0.nc",
        [-10.0, -11.0, 20.0],
        [-50.0, -51.0, -50.0],
        ["2013-01-01 15:00", "2013-01-01 15:00", "2013-01-01 15:01"],
        pressure_weight=WEIGHTS,
    )
    monkeypatch.setattr(prepare, "download_l2", lambda product, day, cache: l2)
    row = prepare.prepare_day("GOSAT-1 v10.0", GOSAT_BRAZIL, "2013-01-01", tmp_path / "month")
    assert row["status"] == "prepared"
    assert (row["l2_soundings"], row["soundings"], row["sources"]) == (3, 2, 32)  # 16 levels

    run = tmp_path / "month" / "20130101"
    with xr.open_dataset(run / "obs_GOSAT-BRAZIL_20130101.nc") as f:
        obs = f.load()
    assert obs.sizes["sounding_id"] == 2
    # Stage 1: the met-only surface-pressure run and a main job that writes its own input.
    sp_text = (run / "surface_pressure" / "input.txt").read_text()
    assert "DummySource" in sp_text
    assert "np.float" not in sp_text
    assert "S0002_PS, Pressure (Pa)" in sp_text
    assert "_T, Temperature" not in sp_text  # surface pressure only
    assert "prepare main-input" in (run / "job.sbatch").read_text()
    assert not (run / "input.txt").exists()

    # Stage 2, once the surface-pressure run has written its output.
    p0 = obs["pressure_levels"].values[:, 0].astype(float) * 100
    model = [p0[0] * 1.01, p0[1] * 0.99]
    for name in ["S0001", "S0002"]:
        text = surface_pressure_output(name, ["S0001", "S0002"], model)
        (run / "surface_pressure" / "output" / f"ReleaseMet_{name}_C1.txt").write_text(text)
    plan = prepare.write_main_input(run)
    groups = pd.read_csv(run / "release_groups.csv", dtype=str)
    assert groups["release_group"].nunique() == 2  # different grid cells: no shared releases (D20)
    surface = pd.read_csv(run / "surface_pressure.csv")
    assert surface["ratio"].round(6).tolist() == [1.01, 0.99]
    # The level-1 layer now starts at the model surface (D19).
    level_1 = plan[plan["level_number"] == 1]
    bottom = level_1["z_pa"] + level_1["dz_pa"] / 2
    assert bottom.to_numpy() == pytest.approx(model)
    saved = pd.read_csv(run / "release_plan.csv")
    assert sorted(saved["sounding_id"].astype(str).unique()) == sorted(obs["sounding_id"].values)
    text = (run / "input.txt").read_text()
    assert text.count("\nSourceID_") == 32
    assert "ReleaseMet_S0002" in text  # release-point met requested in the run (X21)

    # An implausible model surface pressure stops the job.
    sp_file = run / "surface_pressure" / "output" / "ReleaseMet_S0001_C1.txt"
    sp_file.write_text(surface_pressure_output("S0001", ["S0001", "S0002"], [p0[0] * 1.5, p0[1]]))
    with pytest.raises(ValueError, match="Implausible"):
        prepare.write_main_input(run)
    # So does a modelled layer reaching above the domain top (here at 70 hPa; layer 16 tops at
    # about 65 hPa): D21.
    sp_file.write_text(surface_pressure_output("S0001", ["S0001", "S0002"], model, top_pa=7000.0))
    with pytest.raises(ValueError, match="above the domain top"):
        prepare.write_main_input(run)


def test_prepare_day_without_l2_file(tmp_path, monkeypatch):
    def missing(product, day, cache):
        raise L2FileMissing("none")

    monkeypatch.setattr(prepare, "download_l2", missing)
    row = prepare.prepare_day("GOSAT-1 v10.0", GOSAT_BRAZIL, "2013-01-06", tmp_path)
    assert row["status"] == "no L2 file"
