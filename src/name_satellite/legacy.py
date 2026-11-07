"""Readers for files written by the LEGACY ACRG pipeline.

These are used to check the new code against legacy output (see docs/reference-data.md). They are
not part of the new workflow itself.
"""

import pandas as pd

# Column order of a cuboid source line in a NAME III "Sources:" block, as written by
# NAME_scripts/satellite/BackRuns_satellite.scr. Only the columns we use are named.
_SOURCE_COLUMNS = {
    0: "name",
    1: "start_time",
    2: "stop_time",
    3: "particles",
    12: "x",
    13: "y",
    14: "z",
    15: "dx",
    16: "dy",
    17: "dz",
}


def read_name_input_sources(path):
    """Read the particle sources from a legacy NAME input file.

    Each satellite release is a line starting ``SourceID_<n>``, for example::

        SourceID_1, 01/01/2013 14:59:50, 01/01/2013 14:58:50, 30/sec, INERT-TRACER 1.0 g/s,
        Cuboid, Lat-Long, Pa, , , No, No, -34.618, -10.390, 99528.6, 0.0960, 0.0944, 3667.3, ...

    Returns:
        DataFrame with one row per source, in file order, with columns:
        ``name``; ``start_time`` and ``stop_time`` (datetime; for a backwards run start is LATER
        than stop); ``particles_per_second`` (int); ``x``, ``y``, ``dx``, ``dy`` (degrees);
        ``z``, ``dz`` (Pa, the layer centre and thickness).
    """
    rows = []
    with open(path) as f:
        for line in f:
            if not line.startswith("SourceID_"):
                continue
            fields = [field.strip() for field in line.split(",")]
            rows.append({name: fields[i] for i, name in _SOURCE_COLUMNS.items()})
    if not rows:
        raise ValueError(f"No SourceID_ lines found in {path}")

    sources = pd.DataFrame(rows)
    for column in ("x", "y", "z", "dx", "dy", "dz"):
        sources[column] = sources[column].astype(float)
    for column in ("start_time", "stop_time"):
        sources[column] = pd.to_datetime(sources[column], format="%d/%m/%Y %H:%M:%S")
    if not sources["particles"].str.endswith("/sec").all():
        raise ValueError("Expected particle rates of the form '<n>/sec'")
    sources["particles_per_second"] = sources["particles"].str.removesuffix("/sec").astype(int)
    return sources.drop(columns="particles")
