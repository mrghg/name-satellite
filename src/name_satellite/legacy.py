"""Readers for files written by the LEGACY ACRG pipeline.

These are used to check the new code against legacy output (see docs/reference-data.md). They are
not part of the new workflow itself.
"""

import gzip
from pathlib import Path

import numpy as np
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


# Quantity names in legacy met-extraction files -> our column names (with units).
_MET_QUANTITIES = {
    "Temperature (C)": "temperature_c",
    "Boundary layer depth": "boundary_layer_depth_m",
    "Pressure (Pa)": "pressure_pa",
    "Wind speed": "wind_speed_m_s",
    "Wind direction (degrees)": "wind_direction_deg",
}


def read_legacy_release_met(path):
    """Read one legacy release-point met file (from the separate legacy "met" NAME run).

    The legacy run script (``BackRuns_satellite.scr``) runs NAME a second time to sample the met
    at each release point and time. It writes one small file per release, named
    ``Met_<lon>_<lat>_<ID>_C1_X1_Y1.txt[.gz]`` where ``<ID>`` is ``PPPLL``: a 3-digit point
    (sounding) number and a 2-digit level number, both 1-based.

    Returns:
        dict with ``point_number``, ``level_number``, ``time`` (UTC, minute precision), ``lon``,
        ``lat``, ``z_pa`` (the release height requested, Pa), and the met values
        ``temperature_c`` (°C), ``boundary_layer_depth_m``, ``pressure_pa``, ``wind_speed_m_s``,
        ``wind_direction_deg``.

    Note: the release height is given to NAME in Pa, so ``pressure_pa`` is (to rounding) just the
    requested ``z_pa``, i.e. the centre pressure of the release layer, not an independent value.
    The exception is a requested pressure below the model's ground: then NAME samples at the
    surface and ``pressure_pa`` is lower (≈3% of legacy soundings, docs/deviations.md X9).
    """
    path = Path(path)
    release_id = path.name.split("_")[3]
    if len(release_id) != 5 or not release_id.isdigit():
        raise ValueError(f"{path.name}: expected a 5-digit PPPLL release ID")

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        lines = f.read().splitlines()
    rows = [[cell.strip() for cell in line.split(",")] for line in lines]
    heading = next(i for i, row in enumerate(rows) if row[0] == "T")
    fields_line = next(i for i, line in enumerate(lines) if line.strip() == "Fields:")
    quantities = rows[fields_line + 3][3:8]
    if sorted(quantities) != sorted(_MET_QUANTITIES):
        raise ValueError(f"{path.name}: unexpected quantities {quantities}")
    z_labels = set(rows[heading - 2][3:8])
    if len(z_labels) != 1:
        raise ValueError(f"{path.name}: the five quantities are at different heights {z_labels}")
    data = rows[heading + 1 :]
    if len(data) != 1:
        raise ValueError(f"{path.name}: expected one data row, found {len(data)}")
    values = data[0]

    result = {
        "point_number": int(release_id[:3]),
        "level_number": int(release_id[3:]),
        "time": pd.to_datetime(values[0].removesuffix(" UTC"), format="%d/%m/%Y %H:%M"),
        "lon": float(values[1]),
        "lat": float(values[2]),
        "z_pa": float(z_labels.pop().removeprefix("Z = ").removesuffix(" Pa")),
    }
    for quantity, value in zip(quantities, values[3:8], strict=True):
        result[_MET_QUANTITIES[quantity]] = float(value)
    return result


def read_legacy_met_folder(folder):
    """Read every release-point met file in a legacy ``Met/Met_<run name>`` folder.

    Returns:
        DataFrame with one row per release, sorted by (point_number, level_number); see
        ``read_legacy_release_met`` for the columns.
    """
    files = sorted(Path(folder).glob("Met_*_C1_X1_Y1.txt*"))
    if not files:
        raise ValueError(f"No met files in {folder}")
    met = pd.DataFrame([read_legacy_release_met(f) for f in files])
    if met.duplicated(["point_number", "level_number"]).any():
        raise ValueError(f"{folder}: repeated release IDs")
    return met.sort_values(["point_number", "level_number"], ignore_index=True)


def legacy_grid_coordinates(origin, resolution, count):
    """Grid coordinates built exactly as the legacy processing did (``define_grid``, NAME III).

    Legacy uses ``numpy.arange(origin, origin + count*resolution - 0.01*resolution, resolution)``,
    which differs from ``origin + resolution * i`` by up to ~1e-12 degrees. That is irrelevant
    physically but changes the last bits of grid-cell areas and histogram bin edges, so we use the
    same construction when reproducing legacy files bit for bit.

    Args:
        origin, resolution: from the fields-file header (degrees); count: number of cells.
    """
    origin = float(origin)
    resolution = float(resolution)
    coords = np.arange(origin, origin + count * resolution - 0.01 * resolution, resolution)
    if len(coords) != count:
        raise ValueError(f"Expected {count} coordinates, got {len(coords)}")
    return coords
