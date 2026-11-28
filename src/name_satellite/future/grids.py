"""The nested output grids of the future method, and how they are stored.

Where this sits: used by ``future.setup`` (to write each release group's NAME grids) and
``future.processing`` (to place NAME's output into fixed-size blocks). See
docs/future/scoping.md, "Using the nested grids in an inversion" and "Storing the output".

Three levels, all cut from global grids that share the origin (-180°, -90°):

* fine: 0.1° cells, a box of ±5° around each release group;
* medium: 0.5° cells, a box of ±25° around each release group;
* coarse: 2.5° cells, the whole globe (144 x 72), shared by every release.

2.5 = 5 x 0.5 = 25 x 0.1, so every grid line of a coarser level is also a line of the finer ones.
Each box is **snapped outward to the 2.5° lines**: it covers whole 2.5° cells (and so whole 0.5°
cells), and no coarse cell is ever partly inside a finer box. A ±5° box is then 10° or 12.5° wide
(100 or 125 cells), a ±25° box 50° or 52.5° (100 or 105 cells).

Storage: each sounding's fine and medium footprints are kept as fixed-size blocks (the largest a
snapped box can be: 125 x 125 and 105 x 105 cells), zero-padded, with the integer index of the
block's first cell in the global grid of that level (``offset``) and the box's true size
(``shape``). A file then needs no coordinates per sounding.

Limits: boxes are clipped at the poles; a box crossing the dateline raises an error (not needed
for South America; to be handled when a site needs it).
"""

import math
from dataclasses import dataclass

ORIGIN_LON = -180.0
ORIGIN_LAT = -90.0
SNAP_DEG = 2.5


@dataclass(frozen=True)
class Level:
    """One output-grid level: name, cell size and box half-width (degrees), block size (cells)."""

    name: str
    cell_deg: float
    half_width_deg: float | None  # None: the whole globe
    block_cells: int  # fixed block size (square); for the global level, unused


FINE = Level("fine", 0.1, 5.0, 125)
MEDIUM = Level("medium", 0.5, 25.0, 105)
COARSE = Level("coarse", 2.5, None, 0)
LEVELS = (FINE, MEDIUM, COARSE)
COARSE_NX = 144
COARSE_NY = 72


@dataclass(frozen=True)
class Box:
    """A box's outer edges in degrees: west < east, south < north."""

    west: float
    east: float
    south: float
    north: float


def snapped_box(lon, lat, half_width_deg):
    """The box lon/lat ± half_width, widened outward to the nearest 2.5° lines (degrees).

    Clipped at ±90° latitude.

    Raises:
        ValueError if the box would cross the dateline (±180°).
    """
    west = ORIGIN_LON + math.floor((lon - half_width_deg - ORIGIN_LON) / SNAP_DEG) * SNAP_DEG
    east = ORIGIN_LON + math.ceil((lon + half_width_deg - ORIGIN_LON) / SNAP_DEG) * SNAP_DEG
    south = ORIGIN_LAT + math.floor((lat - half_width_deg - ORIGIN_LAT) / SNAP_DEG) * SNAP_DEG
    north = ORIGIN_LAT + math.ceil((lat + half_width_deg - ORIGIN_LAT) / SNAP_DEG) * SNAP_DEG
    if west < -180.0 or east > 180.0:
        raise ValueError(f"Box around ({lon}, {lat}) crosses the dateline: not supported yet")
    return Box(west, east, max(south, -90.0), min(north, 90.0))


def box_shape(box, cell_deg):
    """(ny, nx): number of cells of size cell_deg in the box."""
    return (round((box.north - box.south) / cell_deg), round((box.east - box.west) / cell_deg))


def box_offset(box, cell_deg):
    """(j0, i0): global index (from the origin, 0-based) of the box's first (south-west) cell."""
    return (round((box.south - ORIGIN_LAT) / cell_deg), round((box.west - ORIGIN_LON) / cell_deg))


def name_grid_line(grid_name, box, cell_deg):
    """The ``Horizontal Grids`` row for NAME: X Min / Y Min are the centre of the first cell."""
    ny, nx = box_shape(box, cell_deg)
    x_min = box.west + cell_deg / 2
    y_min = box.south + cell_deg / 2
    return f"{grid_name}, Lat-Long, {nx}, {ny}, {cell_deg}, {cell_deg}, {x_min:.4f}, {y_min:.4f}"


def global_box():
    return Box(-180.0, 180.0, -90.0, 90.0)


def global_centres(cell_deg):
    """(lat, lon) cell centres of the global grid of this cell size (degrees)."""
    ny, nx = round(180 / cell_deg), round(360 / cell_deg)
    lat = [ORIGIN_LAT + cell_deg * (j + 0.5) for j in range(ny)]
    lon = [ORIGIN_LON + cell_deg * (i + 0.5) for i in range(nx)]
    return lat, lon


def offset_from_centres(lat, lon, cell_deg):
    """(j0, i0) of a grid whose first cell centre is (lat[0], lon[0]), checked to lie on the
    global grid of this cell size.

    Raises:
        ValueError if the first cell is not a cell of the global grid (within 1e-3 cell).
    """
    j = (lat[0] - ORIGIN_LAT) / cell_deg - 0.5
    i = (lon[0] - ORIGIN_LON) / cell_deg - 0.5
    if abs(j - round(j)) > 1e-3 or abs(i - round(i)) > 1e-3:
        raise ValueError(f"Grid starting at ({lat[0]}, {lon[0]}) is not on the {cell_deg}° grid")
    return round(j), round(i)
