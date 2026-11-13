"""Check the local reference data are exactly the files documented in docs/reference-data.md."""

import pytest

from name_satellite import __version__
from name_satellite.validation.refdata import REFERENCE_FILES, sha256sum


def test_package_imports():
    assert __version__


def test_reference_keys_unique():
    keys = [r.key for r in REFERENCE_FILES]
    assert len(keys) == len(set(keys))


@pytest.mark.refdata
@pytest.mark.parametrize("ref", REFERENCE_FILES, ids=lambda r: r.key)
def test_reference_checksum(ref, refdata):
    assert sha256sum(refdata(ref.key)) == ref.sha256


@pytest.mark.refdata
@pytest.mark.slow
def test_conc2018_month_manifest(refdata):
    """Every non-met file of the Jan 2013 conc run matches dev/manifests (sha256)."""
    from pathlib import Path

    from name_satellite.validation.refdata import CONC_2018_MANIFEST

    root = refdata("conc2018_input").parents[2]  # .../gosat_brazil_201301
    manifest = Path(__file__).resolve().parents[1] / CONC_2018_MANIFEST
    lines = manifest.read_text().splitlines()
    assert len(lines) == 435
    for line in lines:
        digest, rel = line.split(maxsplit=1)
        assert sha256sum(root / rel) == digest, rel


@pytest.mark.refdata
def test_conc2018_month_met_counts(refdata):
    """Each day has one met file per NAME source (in Met/), for all 31 days."""
    run = refdata("conc2018_input").parents[1]
    for day in range(1, 32):
        d = f"201301{day:02d}"
        tag = f"SOUTHAMERICA_GOSAT-BRAZIL_column_{d}"
        n_src = sum(
            line.startswith("SourceID_")
            for line in (run / "Input_files" / f"BackRun_{tag}.txt").read_text().splitlines()
        )
        n_met = len(list((run / "Met" / f"Met_{tag}").iterdir()))
        assert n_met == n_src, d


@pytest.mark.refdata
@pytest.mark.slow
def test_mr2019_month_manifest(refdata):
    """Every non-met file of the Jan 2013 2019 MR run matches dev/manifests (sha256)."""
    from pathlib import Path

    from name_satellite.validation.refdata import MR_2019_MONTH_MANIFEST

    root = refdata("mr2019_month_fields_20130131").parents[2]  # .../gosat_brazil_201301_mr_2019
    manifest = Path(__file__).resolve().parents[1] / MR_2019_MONTH_MANIFEST
    lines = manifest.read_text().splitlines()
    assert len(lines) == 281
    for line in lines:
        digest, rel = line.split(maxsplit=1)
        assert sha256sum(root / rel) == digest, rel


@pytest.mark.refdata
def test_mr2019_month_met_counts(refdata):
    """Each day of the 2019 MR run has one met file per NAME source, for all 31 days."""
    run = refdata("mr2019_month_fields_20130131").parents[1]
    for day in range(1, 32):
        tag = f"SOUTHAMERICA_GOSAT-BRAZIL_column_201301{day:02d}"
        n_src = sum(
            line.startswith("SourceID_")
            for line in (run / "Input_files" / f"BackRun_{tag}.txt").read_text().splitlines()
        )
        n_met = len(list((run / "Met" / f"Met_{tag}").iterdir()))
        assert n_met == n_src, tag
