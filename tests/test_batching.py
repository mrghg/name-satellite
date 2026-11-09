"""Tests for name_satellite.batching (several days in one NAME run)."""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from name_satellite.batching import combine_release_plans, split_batched_output


def plan(n):
    return pd.DataFrame(
        {"source_name": [f"SourceID_{i}" for i in range(1, n + 1)], "z_pa": np.arange(n) * 1.0}
    )


def test_combine_numbers_sources_over_the_run():
    combined = combine_release_plans({"d1": plan(2), "d2": plan(3)})
    assert combined["source_name"].tolist() == [f"SourceID_{i}" for i in range(1, 6)]
    assert combined["day"].tolist() == ["d1", "d1", "d2", "d2", "d2"]
    assert combined["day_source_name"].tolist() == [
        "SourceID_1", "SourceID_2", "SourceID_1", "SourceID_2", "SourceID_3",
    ]  # fmt: skip
    assert combined["z_pa"].tolist() == [0.0, 1.0, 0.0, 1.0, 2.0]


def test_combine_rejects_misnamed_plans_and_too_many_sources():
    bad = plan(2).iloc[::-1]
    with pytest.raises(ValueError, match="not SourceID_1"):
        combine_release_plans({"d1": bad})
    with pytest.raises(ValueError, match="at most 3000"):
        combine_release_plans({"d1": plan(2000), "d2": plan(1001)})


def test_split_returns_each_day_as_a_one_day_run():
    combined = combine_release_plans({"d1": plan(2), "d2": plan(3)})
    fields = xr.DataArray(
        np.arange(5.0)[:, None, None] * np.ones((5, 2, 2)),
        dims=("source", "lat", "lon"),
        coords={"source": combined["source_name"].to_numpy()},
    )
    particles = pd.DataFrame({"source_number": [1, 3, 3, 5, 2], "lon": [10.0, 30, 31, 50, 20]})

    f2, p2 = split_batched_output(fields, particles, combined, "d2")
    assert f2["source"].values.tolist() == ["SourceID_1", "SourceID_2", "SourceID_3"]
    assert f2.values[:, 0, 0].tolist() == [2.0, 3.0, 4.0]  # global sources 3, 4, 5
    assert p2["source_number"].tolist() == [1, 1, 3]
    assert p2["lon"].tolist() == [30.0, 31.0, 50.0]

    f1, p1 = split_batched_output(fields, particles, combined, "d1")
    assert f1.values[:, 0, 0].tolist() == [0.0, 1.0]
    assert p1["source_number"].tolist() == [1, 2]
