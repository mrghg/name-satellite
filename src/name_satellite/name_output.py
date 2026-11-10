"""Step 5 of the workflow: read the output files written by a NAME run.

Where this sits (see docs/workflow.md)
--------------------------------------
    ... -> 4. run NAME  ->  **5. read NAME output**  ->  6. combine into column footprints

A NAME back-run for one day writes:

* a **fields file**: for every release (NAME "source"), the time-integrated mixing ratio (or air
  concentration) in a near-surface layer on the output grid. Turned into a footprint later.
  Read with ``read_fields``.
* a **particle file**: one line per particle, giving where and when it stopped (left the domain,
  or reached the end of the run). Used for the boundary sensitivities. Read with
  ``read_acrg_particle_file`` (the ACRG-modified NAME format used by the legacy runs).

These readers only *read* and *check*: they do no unit conversion or other processing. Every
release is identified by its NAME source name (e.g. ``SourceID_17``) or source number, never by
the position of a column or line.
"""

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# ---------------------------------------------------------------------------------------------
# Fields files
# ---------------------------------------------------------------------------------------------

# NAME III text fields files (Output Format "IA") have one header row per piece of column
# metadata, between the "Fields:" line and the "X Index, Y Index, ..." line. Checked for the
# NAME 7.2 legacy files (docs/reference-data.md). Empty names mark rows that are always blank.
FIELD_HEADER_ROWS = (
    "species_category",  # e.g. INERT
    "request",  # output request name, e.g. Req_H1Z1_17
    "quantity",  # e.g. Mixing Ratio, Air concentration
    "species",  # e.g. INERT-TRACER
    "units",  # e.g. ppm s, g s / m^3
    "source",  # NAME source name, e.g. SourceID_17  (our release ID)
    "ensemble_averaging",
    "time_averaging",  # e.g. "30day 1hr 0min integral"
    "horizontal_averaging",
    "vertical_averaging",  # e.g. "No vertical averaging" / "Boundary layer average"
    "",
    "",
    "",
    "field_time",  # end of the integration period, e.g. "02/12/2012 14:58 UTC"
    "z_level",  # e.g. "Z = 20.00000 m agl" or "Boundary layer average"
    "",
)
N_PRELIMINARY_COLUMNS = 4  # X Index, Y Index, X (Lat-Long), Y (Lat-Long)


def _open_text(path):
    """Open a plain or gzip-compressed text file."""
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def _read_field_header(path):
    """Read the header of a NAME fields file.

    Returns:
        (info, column_metadata, n_header_lines):
          * ``info``: dict of the "Key: value" lines at the top (run name, grid origin, ...);
          * ``column_metadata``: DataFrame, one row per field column, columns named as in
            FIELD_HEADER_ROWS;
          * ``n_header_lines``: number of lines before the first data row.
    """
    with _open_text(path) as f:
        header_lines = []
        for line in f:
            header_lines.append(line)
            if line.lstrip().startswith("X Index"):
                break
        else:
            raise ValueError(f"No 'X Index' column-heading line in {path}")

    fields_line = next(
        (i for i, line in enumerate(header_lines) if line.strip() == "Fields:"), None
    )
    if fields_line is None:
        raise ValueError(f"No 'Fields:' line in {path}")

    # Lines before "Fields:" are "Key: value" pairs (run name, grid origin, ...).
    info = {}
    for line in header_lines[:fields_line]:
        if ":" in line:
            key, value = line.split(":", 1)
            info[key.strip()] = value.strip()

    # Lines between "Fields:" and "X Index" hold the per-column metadata.
    metadata_rows = [
        [cell.strip() for cell in line.rstrip("\n").split(",")]
        for line in header_lines[fields_line + 1 : -1]
    ]
    if len(metadata_rows) != len(FIELD_HEADER_ROWS):
        raise ValueError(
            f"{path}: expected {len(FIELD_HEADER_ROWS)} column-metadata rows, "
            f"found {len(metadata_rows)}; the NAME header layout has changed"
        )
    n_fields = int(info["Number of field cols"])
    first = N_PRELIMINARY_COLUMNS
    column_metadata = pd.DataFrame(
        {
            name: row[first : first + n_fields]
            for name, row in zip(FIELD_HEADER_ROWS, metadata_rows, strict=True)
            if name
        }
    )
    if len(column_metadata) != n_fields:
        raise ValueError(f"{path}: header has fewer than {n_fields} field columns")
    return info, column_metadata, len(header_lines)


def read_fields(path):
    """Read a NAME III text fields file (time-integrated, Output Format "IA") into memory.

    NAME omits grid rows in which every field is zero; they are restored here as zeros, so the
    result always covers the whole output grid.

    Args:
        path: fields file, plain text or ``.gz``.

    Returns:
        xarray.DataArray ``fields`` with dims ``(source, lat, lon)``:
          * values exactly as written by NAME, in the units given by the ``units`` coordinate
            (e.g. "ppm s" for mixing ratio, "g s / m^3" for air concentration);
          * ``source``: NAME source name of each column (our release ID);
          * per-source coordinates ``request``, ``quantity``, ``units``, ``species``,
            ``time_averaging``, ``vertical_averaging``, ``field_time``, ``z_level``;
          * ``lon``/``lat``: the grid positions as NAME prints them (origin + index × resolution);
          * attrs: the "Key: value" lines from the top of the file (run name, grid, ...).

    Raises:
        ValueError if the header is not as expected, a source name is repeated, a grid index is
        out of range or repeated, or the printed coordinates disagree with the grid definition.
    """
    info, meta, n_header_lines = _read_field_header(path)
    if meta["source"].duplicated().any():
        raise ValueError(f"{path}: repeated source names in the header")

    nx = int(info["X grid size"])
    ny = int(info["Y grid size"])
    lon = float(info["X grid origin"]) + float(info["X grid resolution"]) * np.arange(nx)
    lat = float(info["Y grid origin"]) + float(info["Y grid resolution"]) * np.arange(ny)

    table = pd.read_csv(path, skiprows=n_header_lines, header=None, skipinitialspace=True)
    # Lines end with a trailing comma, which gives one extra empty column.
    n_columns = N_PRELIMINARY_COLUMNS + len(meta)
    if table.shape[1] == n_columns + 1 and table.iloc[:, -1].isna().all():
        table = table.iloc[:, :-1]
    if table.shape[1] != n_columns:
        raise ValueError(f"{path}: expected {n_columns} columns, found {table.shape[1]}")

    ix = table.iloc[:, 0].to_numpy(dtype=int) - 1  # NAME indices are 1-based
    iy = table.iloc[:, 1].to_numpy(dtype=int) - 1
    if ix.min() < 0 or ix.max() >= nx or iy.min() < 0 or iy.max() >= ny:
        raise ValueError(f"{path}: grid index out of range")
    if pd.Series(ix * ny + iy).duplicated().any():
        raise ValueError(f"{path}: a grid cell appears more than once")
    # NAME prints coordinates to 5 decimal places.
    if not (
        np.allclose(table.iloc[:, 2], lon[ix], atol=1e-4)
        and np.allclose(table.iloc[:, 3], lat[iy], atol=1e-4)
    ):
        raise ValueError(f"{path}: printed coordinates disagree with the grid definition")

    values = np.zeros((len(meta), ny, nx))
    values[:, iy, ix] = table.iloc[:, N_PRELIMINARY_COLUMNS:].to_numpy(dtype=float).T

    per_source = {
        name: ("source", meta[name].to_numpy()) for name in meta.columns if name != "source"
    }
    return xr.DataArray(
        values,
        dims=("source", "lat", "lon"),
        coords={"source": meta["source"].to_numpy(), "lat": lat, "lon": lon, **per_source},
        name="fields",
        attrs=info,
    )


# ---------------------------------------------------------------------------------------------
# Particle files
# ---------------------------------------------------------------------------------------------


def read_fields_folder(folder):
    """Read every footprint fields file of one NAME run and combine them by source.

    A run with more releases than fit in one fields file writes several output groups
    (``name_input.fields_output_group``): ``Fields_grid1_C1_T1_*.txt``,
    ``Fields_grid1_2_C1_T1_*.txt``, ... Each is read with ``read_fields``; the result is ordered
    by source number (SourceID_<n>).

    Raises:
        ValueError if no fields file is found, a source appears in more than one file, or the
        files' grids differ.
    """
    files = sorted(Path(folder).glob("Fields_grid1*_C1_T1_*.txt*"))
    if not files:
        raise ValueError(f"No Fields_grid1 files in {folder}")
    parts = [read_fields(path) for path in files]
    for part in parts[1:]:
        if not (
            np.array_equal(part["lat"], parts[0]["lat"])
            and np.array_equal(part["lon"], parts[0]["lon"])
        ):
            raise ValueError(f"Fields files in {folder} have different grids")
    fields = xr.concat(parts, dim="source")
    names = fields["source"].values.tolist()
    if len(set(names)) != len(names):
        raise ValueError(f"A source appears in more than one fields file in {folder}")
    order = np.argsort([int(name.removeprefix("SourceID_")) for name in names])
    fields = fields.isel(source=order)
    fields.attrs = parts[0].attrs
    return fields


def read_acrg_particle_file(path):
    """Read a particle file written by the ACRG-modified NAME (legacy runs).

    The ACRG modification writes one line for every particle when it stops, i.e. when it leaves
    the domain or reaches the end of the run::

          Long    Lat    Ht  Age(hr)   Id
         -59.62 -14.93 19924    0.00  238

    Args:
        path: particle file, plain text or ``.gz``.

    Returns:
        DataFrame with one row per particle:
          * ``lon``, ``lat``: final position (degrees, 2 decimal places);
          * ``height_m_agl``: final height, metres above ground (verified, see
            dev/notes/legacy-processing.md);
          * ``age_hours``: travel time since release (hours, 2 decimal places);
          * ``source_number``: NAME source number, 1-based, i.e. ``SourceID_<source_number>``.

    Raises:
        ValueError if any line does not have exactly five values (the fixed-width format can run
        values together, e.g. longitudes below -99.99).
    """
    with _open_text(path) as f:
        header = f.readline().split()
        if header != ["Long", "Lat", "Ht", "Age(hr)", "Id"]:
            raise ValueError(f"{path}: unexpected header {header}")
        rows = [line.split() for line in f]
    bad = [i for i, row in enumerate(rows) if len(row) != 5]
    if bad:
        raise ValueError(f"{path}: {len(bad)} lines without 5 values (first: line {bad[0] + 2})")

    particles = pd.DataFrame(
        rows, columns=["lon", "lat", "height_m_agl", "age_hours", "source_number"]
    )
    return particles.astype(
        {
            "lon": float,
            "lat": float,
            "height_m_agl": float,
            "age_hours": float,
            "source_number": int,
        }
    )


def particle_counts(particles):
    """Number of particles recorded for each NAME source number (Series indexed by number)."""
    return particles["source_number"].value_counts().sort_index()


def read_name_fate_file(path):
    """Read a standard NAME (8.x) particle-fate file ("Fate Info? = Yes" particle output).

    NAME writes one line per particle when it stops: when it leaves the domain, or when it reaches
    the maximum travel time. After a header ending
    with the column-heading line, each line is::

        Source,Travel Time (hr),X (Lat-Long),Y (Lat-Long),Z (m agl)
        1,3,-33.83381,-10.57067,177.79

    Unlike the legacy ACRG file, the vertical coordinate is whatever the input requested (the
    ``Z-Coord`` of the particle output block), and is named in the heading. With
    ``RunSettings.fate_met`` NAME also writes met columns (U Ambient, ..., Pressure, ...); of these
    only the pressure is kept.

    Args:
        path: fate file, plain text or ``.gz``.

    Returns:
        (particles, z_units): DataFrame with ``source_number`` (1-based NAME source number),
        ``age_hours``, ``lon``, ``lat``, ``z`` (and ``pressure_pa`` if the file has met columns);
        and the vertical coordinate name, e.g. "m agl" or "Pa".
    """
    with _open_text(path) as f:
        lines = f.read().splitlines()
    heading = next(
        (i for i, line in enumerate(lines) if line.strip().startswith("Source,Travel Time")), None
    )
    if heading is None:
        raise ValueError(f"{path}: no 'Source,Travel Time ...' heading line")
    # With "Met? = Yes" NAME adds met columns after Z and ends each line with a comma.
    columns = [c.strip() for c in lines[heading].rstrip().rstrip(",").split(",")]
    expected = ["Source", "Travel Time (hr)", "X (Lat-Long)", "Y (Lat-Long)"]
    if columns[:4] != expected or not columns[4].startswith("Z ("):
        raise ValueError(f"{path}: unexpected columns {columns}")
    if len(columns) != 5 and "Pressure" not in columns:
        raise ValueError(f"{path}: unexpected columns {columns}")
    z_units = columns[4][3:-1]

    rows = [line.rstrip().rstrip(",").split(",") for line in lines[heading + 1 :] if line.strip()]
    bad = [i for i, row in enumerate(rows) if len(row) != len(columns)]
    if bad:
        raise ValueError(f"{path}: {len(bad)} lines without {len(columns)} values")
    table = pd.DataFrame(rows, columns=columns)
    particles = pd.DataFrame(
        {
            "source_number": table["Source"].astype(int),
            "age_hours": table["Travel Time (hr)"].astype(float),
            "lon": table["X (Lat-Long)"].astype(float),
            "lat": table["Y (Lat-Long)"].astype(float),
            "z": table[columns[4]].astype(float),
        }
    )
    if "Pressure" in columns:
        particles["pressure_pa"] = table["Pressure"].astype(float)
    return particles, z_units
