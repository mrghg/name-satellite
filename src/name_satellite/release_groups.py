"""Step 2a of the workflow: group soundings that can share one set of NAME releases (D20).

Where this sits (see docs/workflow.md)
--------------------------------------
    1. soundings (+ 1b. model surface pressure)  ->  **2a. release groups**  ->  2. plan releases
    (one set of releases per group)  ->  ...  ->  6. a column footprint for EVERY sounding, from
    its group's per-level footprints and its OWN pressure weights and averaging kernel

Why
---
GOSAT takes repeated exposures of the same spot seconds apart: in January 2013 (GOSAT-1 v10,
GOSAT-BRAZIL) most soundings come in threes in the same NAME grid cell and minute. Their footprints
are practically identical, so NAME need only be run once for them: about half the releases. Unlike
the legacy binning (``observations.legacy_bin``, X23), nothing is averaged in the observations:
every sounding keeps its own XCH4, uncertainty, averaging kernel, pressure weights and surface
pressure, and the inversion decides how to treat the soundings of a group (e.g. correlated
errors).

Rules
-----
* Soundings are grouped if they fall in the same cell and within ``max_seconds`` (default 60 s)
  of the group's first sounding. Times in GOSAT-1 v10 and GOSAT-2 files are whole minutes rounded
  down (X22), so exposures a few seconds apart can be up to 60 s apart in the files.
* Cells (``Site.release_cell``): "grid" = the NAME output grid's own cells (default); or
  (dlat, dlon) degrees, cells of that size anchored at the output grid's corner; or None = no
  groups (every sounding has its own releases). Legacy bins had the grid's spacing but were
  anchored at the site box's corner, offset from the grid cells.
* A group's release point: mean position and time of its members; mean of the members' release
  levels (after D19 scaling, each with its own model surface pressure); release box = the
  footprint plus the spread of the member positions (as legacy bins); ID = the sounding ID of the
  earliest member (ties: smallest ID), as for legacy bins (D18).

Output
------
No files written here. ``prepare.write_main_input`` writes the membership as
``<day folder>/release_groups.csv`` (one row per sounding: ``sounding_id``, ``release_group``),
which ``processing`` reads to give each sounding its group's per-level results; the group points
go into ``release_plan.csv``, whose ``release_group`` column is the group's ID. A group of one
sounding has that sounding's own ID.
"""

import numpy as np
import pandas as pd
import xarray as xr

from name_satellite.releases import footprint_box_degrees


def cell_size_and_origin(release_cell, grid):
    """Cell size (dlat, dlon) and corner (lat, lon) of the cells, in degrees.

    Args:
        release_cell: "grid", or (dlat, dlon) in degrees.
        grid: the NAME output grid (``name_input.OutputGrid``; x_min/y_min are the centre of the
            first cell).
    """
    corner = (grid.y_min - grid.dy / 2, grid.x_min - grid.dx / 2)
    if release_cell == "grid":
        return (grid.dy, grid.dx), corner
    dlat, dlon = release_cell
    if not (dlat > 0 and dlon > 0):
        raise ValueError(f"Cell size must be positive; got {release_cell}")
    return (float(dlat), float(dlon)), corner


def group_soundings(obs, release_cell, grid, max_seconds=60):
    """Group number (0, 1, ...) of each sounding of ``obs``, in ``obs`` order.

    Args:
        obs: observation Dataset (``observations.to_obs``: ``sounding_id``, ``time``, ``lat``,
            ``lon``).
        release_cell: "grid", (dlat, dlon) or None (each sounding its own group).
        grid: the NAME output grid.
        max_seconds: largest time from a group's first sounding to any other member (s).

    Returns:
        Integer array, one group number per sounding; groups are numbered in order of their
        first sounding's time.

    Raises:
        ValueError if a (dlat, dlon) cell size is not positive.
    """
    n = obs.sizes["sounding_id"]
    times = pd.to_datetime(obs["time"].values)
    if release_cell is None:
        order = np.argsort(times.values, kind="stable")
        groups = np.empty(n, dtype=int)
        groups[order] = np.arange(n)
        return groups

    (dlat, dlon), (lat0, lon0) = cell_size_and_origin(release_cell, grid)
    iy = np.floor((obs["lat"].values - lat0) / dlat).astype(int)
    ix = np.floor((obs["lon"].values - lon0) / dlon).astype(int)
    seconds = (times - times.min()).total_seconds().to_numpy()

    # Walk the soundings in time order; each cell keeps its open group and that group's start time.
    groups = np.empty(n, dtype=int)
    open_group = {}  # (iy, ix) -> (group number, time of its first sounding, s)
    next_group = 0
    for i in np.argsort(seconds, kind="stable"):
        cell = (iy[i], ix[i])
        if cell in open_group and seconds[i] - open_group[cell][1] <= max_seconds:
            groups[i] = open_group[cell][0]
        else:
            open_group[cell] = (next_group, seconds[i])
            groups[i] = next_group
            next_group += 1
    return groups


def release_points(obs, groups, footprint_km, pressure_levels_hpa=None):
    """One release point per group, and the membership of every sounding.

    Args:
        obs: observation Dataset (one entry per sounding).
        groups: from ``group_soundings``.
        footprint_km: sounding footprint diameter (km), for the release box.
        pressure_levels_hpa: optional (sounding, level) array of the levels to release at (e.g.
            scaled to the model surface, D19); default the retrieval's ``pressure_levels``.

    The group's box is the footprint box at the members' mean latitude
    (``footprint_box_degrees``), or, if ``obs`` has ``dlon``/``dlat`` (legacy-binned soundings),
    the largest member box; plus the spread (max - min) of the member positions.

    Returns:
        (points, members):
          * ``points``: Dataset with dim ``release_group`` (the group's ID), in group order, with
            ``time``, ``lon``, ``lat``, ``pressure_levels`` (hPa), ``dlon``, ``dlat`` (release box,
            degrees, full widths) and ``n_soundings``. ``releases.build_release_plan`` takes it.
          * ``members``: DataFrame with ``sounding_id`` and ``release_group`` (the group's ID), one
            row per sounding, in ``obs`` order.
    """
    ids = np.array([str(i) for i in obs["sounding_id"].values])
    times = pd.to_datetime(obs["time"].values)
    lat = obs["lat"].values.astype(float)
    lon = obs["lon"].values.astype(float)
    levels = obs["pressure_levels"].values if pressure_levels_hpa is None else pressure_levels_hpa
    levels = np.asarray(levels, dtype=float)

    rows = []
    group_id = np.empty(len(ids), dtype=object)
    for g in range(groups.max() + 1):
        members = np.flatnonzero(groups == g)
        earliest = members[np.lexsort((ids[members], times[members]))[0]]
        group_id[members] = ids[earliest]
        if "dlat" in obs:
            dlon_fp = obs["dlon"].values[members].max()
            dlat_fp = obs["dlat"].values[members].max()
        else:
            dlon_fp, dlat_fp = footprint_box_degrees(lat[members].mean(), footprint_km)
        rows.append(
            {
                "release_group": ids[earliest],
                "time": times[members].mean(),
                "lon": lon[members].mean(),
                "lat": lat[members].mean(),
                "pressure_levels": levels[members].mean(axis=0),
                "dlon": float(dlon_fp) + (lon[members].max() - lon[members].min()),
                "dlat": float(dlat_fp) + (lat[members].max() - lat[members].min()),
                "n_soundings": len(members),
            }
        )
    dims = "release_group"
    points = xr.Dataset(
        {
            "time": (dims, [r["time"] for r in rows]),
            "lon": (dims, [r["lon"] for r in rows]),
            "lat": (dims, [r["lat"] for r in rows]),
            "pressure_levels": ((dims, "lev"), np.array([r["pressure_levels"] for r in rows])),
            "dlon": (dims, [r["dlon"] for r in rows]),
            "dlat": (dims, [r["dlat"] for r in rows]),
            "n_soundings": (dims, [r["n_soundings"] for r in rows]),
        },
        coords={dims: [r["release_group"] for r in rows]},
    )
    members = pd.DataFrame({"sounding_id": ids, "release_group": group_id.astype(str)})
    return points, members
