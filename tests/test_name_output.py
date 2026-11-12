"""Tests for name_satellite.name_output (workflow step 5: read NAME output)."""

import gzip

import numpy as np
import pandas as pd
import pytest

from name_satellite.legacy import read_name_input_sources
from name_satellite.name_output import (
    particle_counts,
    read_acrg_particle_file,
    read_fields,
    read_fields_folder,
)

# ---------------------------------------------------------------------------------------------
# A tiny hand-made fields file: 3 x 2 grid, 2 sources, one all-zero row left out (as NAME does)
# ---------------------------------------------------------------------------------------------

TINY_FIELDS = """\
 NAME III (version 7.2)
 Run name:                   TINY
 X grid origin:                10.00000
 Y grid origin:                50.00000
 X grid size:                         3
 Y grid size:                         2
 X grid resolution:            0.5000000
 Y grid resolution:            0.2500000
 Number of preliminary cols:            4
 Number of field cols:                2

 Fields:
 , , , , INERT, INERT,
 , , , , Req_1, Req_2,
 , , , , Mixing Ratio, Mixing Ratio,
 , , , , INERT-TRACER, INERT-TRACER,
 , , , , ppm s, ppm s,
 , , , , SourceID_1, SourceID_2,
 , , , , No ensemble averaging, No ensemble averaging,
 , , , , 30day integral, 30day integral,
 , , , , No horizontal averaging, No horizontal averaging,
 , , , , No vertical averaging, No vertical averaging,
 , , , , , ,
 , , , , , ,
 , , , , , ,
 , , , , 02/12/2012 14:58 UTC, 02/12/2012 14:58 UTC,
 , , , , Z = 20.00000 m agl, Z = 20.00000 m agl,
 , , , , , ,
 X Index, Y Index, X (Lat-Long), Y (Lat-Long), , ,
 1, 1, 10.00000, 50.00000, 1.0E+00, 0.0E+00,
 3, 1, 11.00000, 50.00000, 2.5E-01, 4.0E+00,
 2, 2, 10.50000, 50.25000, 0.0E+00, 7.0E-03,
"""


@pytest.fixture
def tiny_fields_file(tmp_path):
    path = tmp_path / "Fields_tiny.txt.gz"
    with gzip.open(path, "wt") as f:
        f.write(TINY_FIELDS)
    return path


def test_tiny_fields_values_and_metadata(tiny_fields_file):
    fields = read_fields(tiny_fields_file)
    assert fields.dims == ("source", "lat", "lon")
    assert fields["source"].values.tolist() == ["SourceID_1", "SourceID_2"]
    np.testing.assert_allclose(fields["lon"], [10.0, 10.5, 11.0])
    np.testing.assert_allclose(fields["lat"], [50.0, 50.25])
    # [source, lat, lon]; every cell not listed in the file is zero.
    expected = np.array(
        [
            [[1.0, 0.0, 0.25], [0.0, 0.0, 0.0]],
            [[0.0, 0.0, 4.0], [0.0, 7e-3, 0.0]],
        ]
    )
    np.testing.assert_array_equal(fields.values, expected)
    assert fields["units"].values.tolist() == ["ppm s", "ppm s"]
    assert fields["request"].values.tolist() == ["Req_1", "Req_2"]
    assert fields.attrs["Run name"] == "TINY"


def test_fields_folder_combines_output_groups_by_source(tmp_path):
    # A run split over two output groups (name_input.fields_output_group): sources 1-2 and 3-4.
    (tmp_path / "Fields_grid1_C1_T1_201212021458.txt").write_text(TINY_FIELDS)
    second = TINY_FIELDS.replace("SourceID_1, SourceID_2", "SourceID_3, SourceID_4")
    (tmp_path / "Fields_grid1_2_C1_T1_201212021458.txt").write_text(second)
    fields = read_fields_folder(tmp_path)
    assert fields["source"].values.tolist() == [f"SourceID_{n}" for n in (1, 2, 3, 4)]
    np.testing.assert_array_equal(fields.sel(source="SourceID_3"), fields.sel(source="SourceID_1"))
    # The same source in two files is an error.
    (tmp_path / "Fields_grid1_3_C1_T1_201212021458.txt").write_text(TINY_FIELDS)
    with pytest.raises(ValueError, match="more than one fields file"):
        read_fields_folder(tmp_path)


def test_tiny_fields_wrong_coordinate_is_rejected(tmp_path):
    path = tmp_path / "bad.txt"
    path.write_text(TINY_FIELDS.replace("3, 1, 11.00000", "3, 1, 11.50000"))
    with pytest.raises(ValueError, match="coordinates disagree"):
        read_fields(path)


def test_tiny_fields_repeated_cell_is_rejected(tmp_path):
    path = tmp_path / "bad.txt"
    path.write_text(TINY_FIELDS + " 1, 1, 10.00000, 50.00000, 1.0E+00, 0.0E+00,\n")
    with pytest.raises(ValueError, match="more than once"):
        read_fields(path)


def test_tiny_fields_repeated_source_is_rejected(tmp_path):
    path = tmp_path / "bad.txt"
    path.write_text(TINY_FIELDS.replace("SourceID_1, SourceID_2", "SourceID_1, SourceID_1"))
    with pytest.raises(ValueError, match="repeated source"):
        read_fields(path)


# ---------------------------------------------------------------------------------------------
# Legacy reference files
# ---------------------------------------------------------------------------------------------


@pytest.mark.refdata
def test_mr2019_fields(refdata):
    fields = read_fields(refdata("mr2019_fields"))
    sources = read_name_input_sources(refdata("mr2019_input"))
    assert fields.shape == (425, 357, 190)
    # Columns are identified by source name; they match the input file's sources.
    assert fields["source"].values.tolist() == sources["name"].tolist()
    assert set(fields["quantity"].values) == {"Mixing Ratio"}
    assert set(fields["units"].values) == {"ppm s"}
    assert set(fields["z_level"].values) == {"Z = 20.00000 m agl"}
    assert (fields.values >= 0).all()
    # Every surface-layer release (level 1) influences the surface. Higher levels can have a zero
    # footprint (particles leave the domain without reaching the surface layer); see
    # dev/notes/legacy-processing.md.
    totals = fields.sum(dim=("lat", "lon")).values.reshape(25, 17)
    assert (totals[:, 0] > 0).all()


@pytest.mark.refdata
def test_conc2018_fields_and_boundary_layer_fields(refdata):
    run = refdata("conc2018_input").parents[1]
    fields = read_fields(
        run / "Fields_files" / "Fields_SOUTHAMERICA_GOSAT-BRAZIL_column_20130101.txt.gz"
    )
    bl = read_fields(
        run / "Fields_files_BL" / "Fields_BL_SOUTHAMERICA_GOSAT-BRAZIL_column_20130101.txt.gz"
    )
    assert set(fields["units"].values) == {"g s / m^3"}
    assert set(bl["vertical_averaging"].values) == {"Boundary layer average"}
    assert fields["source"].values.tolist() == bl["source"].values.tolist()


@pytest.mark.refdata
def test_mr2019_particles(refdata):
    particles = read_acrg_particle_file(refdata("mr2019_particles"))
    sources = read_name_input_sources(refdata("mr2019_input"))
    assert len(particles) == 562_925
    # A 60 s release at N particles/s gives 60 * N + 1 particles (NAME releases one at each end).
    expected = pd.Series(
        60 * sources["particles_per_second"].to_numpy() + 1, index=np.arange(1, len(sources) + 1)
    )
    pd.testing.assert_series_equal(particle_counts(particles), expected, check_names=False)
    # 721 h run: no particle is older than the run.
    assert particles["age_hours"].between(0, 721).all()


@pytest.mark.refdata
def test_conc2018_reprocessed_particle_file_is_the_same_data(refdata):
    original = read_acrg_particle_file(refdata("conc2018_particles"))
    reprocessed_path = (
        refdata("conc2018_particles").parents[1]
        / "Reprocessed"
        / "particle_location_SOUTHAMERICA_GOSAT-BRAZIL_column_20130101.txt.gz"
    )
    pd.testing.assert_frame_equal(original, read_acrg_particle_file(reprocessed_path))


@pytest.mark.refdata
@pytest.mark.slow
def test_conc2018_all_january_fields_and_particles_match_sources(refdata):
    """For every day: fields columns and particle counts match the input file's sources."""
    run = refdata("conc2018_input").parents[1]
    for day in range(1, 32):
        tag = f"SOUTHAMERICA_GOSAT-BRAZIL_column_201301{day:02d}"
        sources = read_name_input_sources(run / "Input_files" / f"BackRun_{tag}.txt")
        fields = read_fields(run / "Fields_files" / f"Fields_{tag}.txt.gz")
        assert fields["source"].values.tolist() == sources["name"].tolist(), tag
        counts = particle_counts(
            read_acrg_particle_file(run / "Particle_files" / f"particle_location_{tag}.txt.gz")
        )
        assert counts.tolist() == (60 * sources["particles_per_second"] + 1).tolist(), tag


# ---------------------------------------------------------------------------------------------
# Standard NAME 8.x particle-fate files
# ---------------------------------------------------------------------------------------------

TINY_FATE = """\
 NAME III (version 8.7.1)
 Run name:            debug_20130101
 Backward trajectories

 Set of particle/puff information: Particles
 Source,Travel Time (hr),X (Lat-Long),Y (Lat-Long),Z (m agl)
 1,3,-33.83381,-10.57067,177.79
 2,2.5,-33.81963,-10.54067,383.6606
 1,0,-34.6,-10.4,20000
"""


def test_tiny_fate_file(tmp_path):
    from name_satellite.name_output import read_name_fate_file

    path = tmp_path / "Particles_C1.txt"
    path.write_text(TINY_FATE)
    particles, z_units = read_name_fate_file(path)
    assert z_units == "m agl"
    assert particles["source_number"].tolist() == [1, 2, 1]
    assert particles["age_hours"].tolist() == [3.0, 2.5, 0.0]
    assert particles["z"].iloc[1] == pytest.approx(383.6606)
    assert particle_counts(particles).to_dict() == {1: 2, 2: 1}


# As NAME 8.7.1 writes it with "Met? = Yes" (RunSettings.fate_met): met columns after Z, and a
# trailing comma on every line (cut down to two met columns plus Pressure, less padding).
TINY_FATE_WITH_MET = """\
 NAME III (version 8.7.1)
 Run name:            fate_test
 Backward trajectories

 Set of particle/puff information: Particles
 Source,Travel Time (hr),X (Lat-Long),Y (Lat-Long),Z (m agl),  U Ambient,  V Ambient,  Pressure,
 1,43.89711,-19.42067,-16.07199,102.2506,       9.662325,  -3.843483,  100489.8,
 2,44.09135,-19.42704,-12.39873,1913.624,       12.67187,  0.3975459,  81052.22,
"""


def test_fate_file_with_met_columns_keeps_pressure(tmp_path):
    from name_satellite.name_output import read_name_fate_file

    path = tmp_path / "Particles_C1.txt"
    path.write_text(TINY_FATE_WITH_MET)
    particles, z_units = read_name_fate_file(path)
    assert z_units == "m agl"
    assert particles["source_number"].tolist() == [1, 2]
    assert particles["pressure_pa"].tolist() == [100489.8, 81052.22]
    assert particles["z"].tolist() == [102.2506, 1913.624]


# As the ACRG-patched NAME writes it (met.PATCHED_NAME_DIR, X20): particle index and fate code
# after Z, here with met columns too.
TINY_FATE_PATCHED = """\
 NAME III (version 8.7.1)
 Run name:            fate_test

 Set of particle/puff information: Particles
 Source,Travel Time (hr),X (Lat-Long),Y (Lat-Long),Z (m agl),Particle Index,Fate Code,  Pressure,
 1,43.89711,-19.42067,-16.07199,102.2506,17,8,  100489.8,
 2,720,-25.1,-12.3,1913.624,1801,16,  81052.22,
"""


def test_fate_file_from_patched_name(tmp_path):
    from name_satellite.name_output import read_name_fate_file

    path = tmp_path / "Particles_C1.txt"
    path.write_text(TINY_FATE_PATCHED)
    particles, _ = read_name_fate_file(path)
    assert particles["particle_index"].tolist() == [17, 1801]
    assert particles["fate_code"].tolist() == [8, 16]
    assert particles["pressure_pa"].tolist() == [100489.8, 81052.22]


def test_drop_repeated_fate_records():
    from name_satellite.name_output import drop_repeated_fate_records

    particles = pd.DataFrame(
        {
            "source_number": [1, 1, 1, 2, 2],
            "particle_index": [7, 8, 7, 7, 7],
            # particle (1, 7) twice, 10 s apart; particle (2, 7) twice, 180 s apart (30-day run)
            "age_hours": [47.99717, 10.0, 47.99995, 51.03133, 50.98133],
            "fate_code": [8, 16, 8, 8, 8],
        }
    )
    kept = drop_repeated_fate_records(particles)
    assert kept.index.tolist() == [0, 1, 4]  # the earlier record of each pair
    particles.loc[2, "age_hours"] = 49.0  # an hour later: longer than a sync step
    with pytest.raises(ValueError, match="X20"):
        drop_repeated_fate_records(particles)
