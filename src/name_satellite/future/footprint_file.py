"""Read a future-method footprint file (one month) and put its footprints on the global grids.

Where this sits (layout: docs/future/file-format.md)
----------------------------------------------------
    ``future.setup`` -> NAME -> ``future.processing`` (writes the month file) -> **this module**
    -> the inversion, or ``notebooks/future_footprints.ipynb`` (plots)

The file stores each sounding's column footprint at three levels: a 0.1° box (fine), a 0.5° box
(medium), both as zero-padded fixed-size blocks with their place in the global grid of that level,
and the global 2.5° grid (coarse). Footprints are in (mol/mol) / (mol m-2 s-1): the column
enhancement from a flux of 1 mol m-2 s-1 in that cell, so summing cells gives the response to a
uniform flux, and a cell's value can be **spread evenly** over the finer cells inside it (each
gets value / n_subcells) without changing any total.

Because every box is snapped to the 2.5° lines, each coarse cell is wholly inside or wholly
outside a finer box, and box offsets and shapes are whole numbers of coarser cells (checked).

What this module gives:

* ``open_footprints(path)``: the dataset, checked for the variables this module needs;
* ``sounding_footprint(ds, sounding_id, level)``: one sounding's footprint on one level, cut to
  its box, with latitude and longitude;
* ``box_edges(ds, sounding_id, level)``: the box's edges in degrees;
* ``finest_footprint(ds, level, sounding_ids)``: the sum of soundings' footprints on the global
  grid of ``level``, taking for each cell the finest level a sounding has there and spreading
  coarser values evenly (how an inversion would use the nested grids);
* ``aggregate(footprint, factor)``: sum blocks of factor x factor cells (to compare levels);
* ``footprint_density(footprint)``: footprint per unit area, which can be compared (and plotted)
  across levels;
* ``sounding_table(ds)``: one row per sounding: place, time, how much of its response lies in its
  boxes, and the share of its particles that left through the top.
"""

from pathlib import Path

import numpy as np
import xarray as xr

from name_satellite.footprints import grid_cell_area_m2
from name_satellite.future.grids import COARSE, FINE, MEDIUM, Box

LEVELS = {"fine": FINE, "medium": MEDIUM, "coarse": COARSE}

REQUIRED_VARIABLES = (
    "srr_fine", "fine_offset", "fine_shape",
    "srr_medium", "medium_offset", "medium_shape",
    "srr_coarse", "end_positions", "column_weight",
    "time", "lat", "lon", "release_group", "release_lat", "release_lon", "release_time",
    "end_time", "top_exit_mean_age",
)  # fmt: skip


def open_footprints(path):
    """Open a month file (lazily) and check it has the future-method layout.

    path: the ``<SITE>_<PRODUCT>_<YYYYMM>_future.nc`` file written by ``future.processing``.
    Returns an xarray.Dataset indexed by ``sounding_id``. Raises if a variable is missing.
    """
    ds = xr.open_dataset(Path(path))
    missing = [name for name in REQUIRED_VARIABLES if name not in ds.variables]
    if missing:
        raise ValueError(f"{path} is not a future-method footprint file: missing {missing}")
    return ds


def _cells_per(coarser, finer):
    """How many ``finer`` cells span one ``coarser`` cell along each axis (2.5/0.5 = 5, etc.)."""
    ratio = LEVELS[coarser].cell_deg / LEVELS[finer].cell_deg
    if abs(ratio - round(ratio)) > 1e-9:
        raise ValueError(f"{coarser} cells are not whole numbers of {finer} cells")
    return round(ratio)


def _box_indices(one, level):
    """Global (lat, lon) index of the box's first cell, and its shape (one sounding's Dataset)."""
    (j0, i0), (ny, nx) = one[f"{level}_offset"].values, one[f"{level}_shape"].values
    return int(j0), int(i0), int(ny), int(nx)


def sounding_footprint(ds, sounding_id, level):
    """One sounding's footprint on one level as a DataArray with ``lat``/``lon`` cell centres.

    level: "fine" (0.1° box), "medium" (0.5° box) or "coarse" (global 2.5°). The fine and medium
    blocks are cut to their box (the zero padding is dropped).
    Returns (mol/mol) / (mol m-2 s-1), dimensions (lat, lon), degrees north / east.
    """
    one = ds.sel(sounding_id=sounding_id)
    if level == "coarse":
        return one["srr_coarse"].rename(coarse_lat="lat", coarse_lon="lon").load()
    j0, i0, ny, nx = _box_indices(one, level)
    block = one[f"srr_{level}"].values[:ny, :nx]
    lat = ds[f"{level}_lat_global"].values[j0 : j0 + ny]
    lon = ds[f"{level}_lon_global"].values[i0 : i0 + nx]
    return xr.DataArray(block, dims=("lat", "lon"), coords={"lat": lat, "lon": lon},
                        attrs={"units": ds["srr_coarse"].attrs.get("units", "")})  # fmt: skip


def box_edges(ds, sounding_id, level):
    """Edges of a sounding's box on ``level`` (degrees) as a ``Box``; the globe for "coarse"."""
    if level == "coarse":
        return Box(-180.0, 180.0, -90.0, 90.0)
    footprint = sounding_footprint(ds, sounding_id, level)
    half = LEVELS[level].cell_deg / 2
    lat, lon = footprint["lat"].values, footprint["lon"].values
    return Box(west=round(lon[0] - half, 6), east=round(lon[-1] + half, 6),
               south=round(lat[0] - half, 6), north=round(lat[-1] + half, 6))  # fmt: skip


def _global_coords(ds, level):
    """Cell-centre latitudes and longitudes of the global grid of ``level`` (degrees)."""
    if level == "coarse":
        return ds["coarse_lat"].values, ds["coarse_lon"].values
    return ds[f"{level}_lat_global"].values, ds[f"{level}_lon_global"].values


def _spread(values, factor):
    """Spread each cell's value evenly over factor x factor finer cells (totals unchanged)."""
    return np.repeat(np.repeat(values, factor, axis=-2), factor, axis=-1) / factor**2


def finest_footprint(ds, level="medium", sounding_ids=None):
    """Sum of soundings' footprints on the global grid of ``level``, finest level first.

    For each sounding and each cell of the global ``level`` grid, the value comes from the finest
    level the sounding has there: its fine box, else its medium box, else the coarse grid; a
    coarser value is spread evenly over the finer cells inside it. This is the "use the finest
    level available" rule of docs/future/scoping.md.

    level: "coarse", "medium" or "fine" (the fine global grid is 1800 x 3600: ~50 MB).
    sounding_ids: the soundings to sum (default: all). One ID gives that sounding's footprint.
    Returns a DataArray (lat, lon), (mol/mol) / (mol m-2 s-1). Its total is close to, but not
    exactly, the sum of the soundings' coarse footprints: finer cells sum to coarser ones only
    approximately (docs/future/file-format.md, "Grids").

    Computed without building a global array per sounding: start from the summed coarse grid
    spread to ``level``, then for each sounding swap the spread coarse values inside its medium box
    for its medium block (and, for "fine", the spread medium values inside its fine box for its
    fine block).
    """
    if sounding_ids is None:
        sounding_ids = ds["sounding_id"].values
    sub = ds.sel(sounding_id=np.atleast_1d(sounding_ids))
    lat, lon = _global_coords(ds, level)
    blocks = {
        name: sub[f"srr_{name}"].values.astype(float) for name in ("coarse", "medium", "fine")
    }
    offsets = {name: sub[f"{name}_offset"].values for name in ("medium", "fine")}
    shapes = {name: sub[f"{name}_shape"].values for name in ("medium", "fine")}

    total = _spread(blocks["coarse"].sum(axis=0), _cells_per("coarse", level))
    finer_levels = {"coarse": [], "medium": ["medium"], "fine": ["medium", "fine"]}[level]
    coarser_of = {"medium": "coarse", "fine": "medium"}
    for box_level in finer_levels:
        coarser = coarser_of[box_level]
        to_level = _cells_per(box_level, level)  # box cells -> cells of the output grid
        per_coarser = _cells_per(coarser, box_level)  # coarser cells -> box cells
        if (offsets[box_level] % per_coarser).any() or (shapes[box_level] % per_coarser).any():
            raise ValueError(f"A {box_level} box is not aligned with the {coarser} cells")
        for k in range(len(blocks["coarse"])):
            (j0, i0), (ny, nx) = offsets[box_level][k], shapes[box_level][k]
            # The box in units of the coarser level's global cells.
            cj, ci = j0 // per_coarser, i0 // per_coarser
            cny, cnx = ny // per_coarser, nx // per_coarser
            if coarser == "medium":  # index into the sounding's medium block, not the globe
                cj, ci = cj - offsets["medium"][k][0], ci - offsets["medium"][k][1]
            old = blocks[coarser][k, cj : cj + cny, ci : ci + cnx]
            new = blocks[box_level][k, :ny, :nx]
            # Swap the coarser values (spread) for the box's own values, at the output resolution.
            change = _spread(new, to_level) - _spread(_spread(old, per_coarser), to_level)
            oj, oi = j0 * to_level, i0 * to_level
            total[oj : oj + ny * to_level, oi : oi + nx * to_level] += change
    return xr.DataArray(total, dims=("lat", "lon"), coords={"lat": lat, "lon": lon},
                        attrs={"units": ds["srr_coarse"].attrs.get("units", "")})  # fmt: skip


def aggregate(footprint, factor):
    """Sum blocks of factor x factor cells of a (lat, lon) footprint (e.g. 0.1° -> 0.5°: 5).

    The footprint's extent must be a whole number of blocks (true for the snapped boxes). The new
    cell centres are the means of the old ones.
    """
    if footprint.sizes["lat"] % factor or footprint.sizes["lon"] % factor:
        raise ValueError(f"Footprint of shape {footprint.shape} is not a whole number of blocks")
    return footprint.coarsen(lat=factor, lon=factor).sum()


def footprint_density(footprint):
    """A (lat, lon) footprint divided by its cells' areas: (mol/mol) / (mol s-1).

    The column enhancement from a point source of 1 mol s-1 in the cell. Unlike the footprint per
    cell, it does not depend on the cell size, so the three levels can be compared and drawn on
    one colour scale. Areas as in the processing (``footprints.grid_cell_area_m2``).
    """
    area = grid_cell_area_m2(footprint["lat"].values, footprint["lon"].values)
    density = footprint / area
    density.attrs = {"units": "(mol/mol) / (mol s-1)"}
    return density


def sounding_table(ds):
    """One row per sounding (index ``sounding_id``): its place and time, and summary numbers.

    Columns: time, lat, lon, release_group, release_lat, release_lon, release_time, end_time,
    top_exit_mean_age (h), column_weight (Σ over modelled levels), total_response (Σ of the coarse
    footprint, (mol/mol)/(mol m-2 s-1): the enhancement from a uniform 1 mol m-2 s-1 flux),
    in_fine_box and in_medium_box (fraction of total_response inside each box, from the coarse
    cells the box covers), left_top (fraction of the column's particles that left through the
    20 km top; the rest reached the maximum travel time below it).
    """
    columns = ["time", "lat", "lon", "release_group", "release_lat", "release_lon",
               "release_time", "end_time", "top_exit_mean_age"]  # fmt: skip
    table = ds[columns].to_dataframe()[columns]
    weight = ds["column_weight"].sum("release_level").values
    coarse = ds["srr_coarse"].values.astype(float)
    total = coarse.sum(axis=(1, 2))
    table["column_weight"] = weight
    table["total_response"] = total
    for level in ("fine", "medium"):
        per_coarse = _cells_per("coarse", level)
        offsets, shapes = ds[f"{level}_offset"].values, ds[f"{level}_shape"].values
        if (offsets % per_coarse).any() or (shapes % per_coarse).any():
            raise ValueError(f"{level} boxes are not aligned with the 2.5° cells")
        inside = np.empty(len(total))
        boxes = zip(offsets // per_coarse, shapes // per_coarse, strict=True)
        for k, ((j0, i0), (ny, nx)) in enumerate(boxes):
            inside[k] = coarse[k, j0 : j0 + ny, i0 : i0 + nx].sum()
        table[f"in_{level}_box"] = inside / total
    top = ds["end_positions"].isel(end_height=-1).sum(("coarse_lat", "coarse_lon")).values
    table["left_top"] = top / weight
    return table
