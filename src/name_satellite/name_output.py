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


def read_acrg_particle_file(path):
    """Read a particle file written by the ACRG-modified NAME (legacy runs).

    The ACRG modification writes one line for every particle when it stops, i.e. when it leaves
    the domain or reaches the end of the run (``Particle.F90``, ``InactivateParticle``)::

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
