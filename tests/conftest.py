import pytest

from name_satellite.validation.refdata import ref_path


@pytest.fixture
def refdata():
    """Return a function mapping a reference-file key to its path, skipping if absent."""

    def _get(key: str):
        path = ref_path(key)
        if not path.exists():
            pytest.skip(f"reference data not available: {path}")
        return path

    return _get
