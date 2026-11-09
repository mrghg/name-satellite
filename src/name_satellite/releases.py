"""Step 2 of the workflow: plan the NAME particle releases for each satellite sounding.

Where this sits (see docs/workflow.md)
--------------------------------------
    1. read soundings  ->  **2. plan releases**  ->  3. write NAME inputs  ->  ...

A satellite retrieval reports its vertical sensitivity on a set of pressure *levels* (for GOSAT,
20 levels from the surface up to 0.1 hPa). NAME cannot release particles at a level; it releases
them uniformly through a *layer* (a box with a top and a bottom). This module turns each
sounding's levels into layers, one per level, and decides how many particles each one releases.

How the layers are built (identical to the legacy ACRG code)
------------------------------------------------------------
Legacy source: ``acrg/satellite/common.py``, ``midpoint_bounds(set_as_edges=False)``.

Levels are numbered from the surface upwards, with pressures p[0] > p[1] > ... > p[n-1], where
p[0] is the surface pressure. The layer boundaries are:

    * bottom of layer 0      = p[0]                       (the surface)
    * boundary between i-1, i = (p[i-1] + p[i]) / 2       (halfway between the two levels)
    * top of the last layer   = p[n-1] - (p[n-2] - p[n-1]) / 2, but never below 0 hPa

So each layer extends halfway to its neighbouring levels. Layer 0 is only half as thick as the
others because it starts at the surface. This is consistent with the retrieval's pressure weights,
which give level 0 roughly half weight.

NAME is told the *centre* and *thickness* of each layer, in Pa.

Legacy option: NAME surface pressure
------------------------------------
The 2013 legacy runs replaced p[0] by a surface pressure taken from a separate NAME run
(``use_name_pressure=True``); see docs/deviations.md X9. Pass ``surface_pressure_hpa`` to
reproduce that. By default the retrieval's own p[0] is used.
"""

import numpy as np
import pandas as pd

HPA_TO_PA = 100.0


def layer_boundaries_hpa(pressure_levels_hpa, surface_pressure_hpa=None):
    """Pressure boundaries (hPa) of the release layers for ONE sounding.

    Args:
        pressure_levels_hpa: retrieval level pressures in hPa, ordered from the surface upwards
            (strictly decreasing). Level 0 is the surface.
        surface_pressure_hpa: optional replacement for level 0 (legacy ``use_name_pressure``).

    Returns:
        Array of n+1 boundaries in hPa, from the surface (index 0) to the top (index n).
        Layer i lies between boundaries i (bottom) and i+1 (top).
    """
    p = np.array(pressure_levels_hpa, dtype=float)  # a copy, so the input is never changed
    if p.ndim != 1 or p.size < 2:
        raise ValueError("Need a 1-D array of at least two pressure levels")
    if surface_pressure_hpa is not None:
        p[0] = surface_pressure_hpa
    if not np.all(np.diff(p) < 0):
        raise ValueError(f"Pressure levels must decrease upwards; got {p}")

    boundaries = np.empty(p.size + 1)
    boundaries[0] = p[0]  # the surface
    boundaries[1:-1] = (p[:-1] + p[1:]) / 2.0  # halfway between neighbouring levels
    top_spacing = p[-2] - p[-1]
    boundaries[-1] = max(p[-1] - top_spacing / 2.0, 0.0)  # extend the top layer symmetrically
    return boundaries


def release_layers_pa(pressure_levels_hpa, surface_pressure_hpa=None):
    """Centre and thickness (both in Pa) of the release layer for each level of ONE sounding.

    These are the ``Z`` and ``dZ`` values of a NAME cuboid source with a ``Pa`` Z-coordinate.

    Args:
        pressure_levels_hpa: retrieval level pressures in hPa, surface first (see
            ``layer_boundaries_hpa``).
        surface_pressure_hpa: optional replacement for level 0 (legacy ``use_name_pressure``).

    Returns:
        (centre_pa, thickness_pa): two arrays with one value per level.
    """
    boundaries = layer_boundaries_hpa(pressure_levels_hpa, surface_pressure_hpa)
    bottom = boundaries[:-1]
    top = boundaries[1:]
    centre_pa = (bottom + top) / 2.0 * HPA_TO_PA
    thickness_pa = (bottom - top) * HPA_TO_PA
    return centre_pa, thickness_pa


def particles_per_second(level_number, level_split=12, low_rate=30, high_rate=3):
    """Particle release rate for a level (legacy rule from ``BackRuns_satellite.scr``).

    Levels at or below ``level_split`` release ``low_rate`` particles per second; higher levels
    release ``high_rate`` (they are assumed to contribute little to surface sensitivity).
    The legacy launcher (``Satellite_launch.scr``) passes ``level_split=12``.

    Args:
        level_number: 1-based level number (1 = the surface layer), as in legacy source IDs.
    """
    if level_number < 1:
        raise ValueError("level_number is 1-based")
    return low_rate if level_number <= level_split else high_rate


# ---------------------------------------------------------------------------------------------
# Horizontal size of a release (the satellite's ground footprint)
# ---------------------------------------------------------------------------------------------

EARTH_RADIUS_KM = 6371.0  # value used by legacy acrg.satellite.common


def footprint_box_degrees(lat, diameter_km=10.5):
    """Width (dlon) and height (dlat), in degrees, of a square box ``diameter_km`` across.

    Legacy (``acrg.satellite.common.distance_lat`` / ``distance_lon``, haversine on a sphere of
    radius 6371 km): GOSAT's ground footprint is ~10.5 km across, and NAME releases particles
    uniformly over a box of this size centred on the sounding. For legacy soundings that were
    averaged into bins, the box was widened by the spread of the averaged positions; that spread is
    not stored in the observation file.

    Returns:
        (dlon, dlat) in degrees; full widths, as NAME's dX and dY.
    """
    half_angle = diameter_km / (2.0 * EARTH_RADIUS_KM)
    a = np.sin(half_angle) ** 2
    dlat = np.degrees(np.abs(2.0 * np.arcsin(np.sqrt(a))))
    dlon = np.degrees(np.abs(2.0 * np.arcsin(np.sqrt(a) / np.cos(np.radians(lat)))))
    return dlon, dlat


# ---------------------------------------------------------------------------------------------
# The release plan: every NAME release of one day, in one table
# ---------------------------------------------------------------------------------------------


def build_release_plan(
    obs,
    n_levels=17,
    release_seconds=60,
    level_split=12,
    footprint_diameter_km=10.5,
    surface_pressure_hpa=None,
    box_degrees=None,
):
    """One row per NAME release (source) for one day of soundings.

    Args:
        obs: observation Dataset in the legacy format (dims ``time``, ``lev``; variables ``lon``,
            ``lat``, ``pressure_levels`` in hPa), one entry per sounding, in time order.
        n_levels: retrieval levels released per sounding, from the surface (legacy 17).
        release_seconds: length of each release (legacy 60 s).
        level_split, footprint_diameter_km: see ``particles_per_second``, ``footprint_box_degrees``.
        surface_pressure_hpa: optional per-sounding replacement of level 0 (legacy X9).
        box_degrees: optional (dlon, dlat) arrays per sounding, overriding the footprint box (e.g.
            the widened boxes of legacy binned soundings).

    Returns:
        DataFrame with columns
          * ``source_name``: "SourceID_<n>", n = 1, 2, ... in sounding order, levels within each;
          * ``point_number``, ``level_number``: 1-based sounding and level numbers;
          * ``release_end``: the sounding time truncated to the whole second (legacy, UTC). In
            NAME's backwards run this is the source's "Stop Time";
          * ``release_start``: ``release_end`` + ``release_seconds`` (NAME's "Start Time");
          * ``lon``, ``lat``, ``dlon``, ``dlat``: box centre and full widths (degrees);
          * ``z_pa``, ``dz_pa``: layer centre and thickness (Pa);
          * ``particles_per_second``.
    """
    times = pd.to_datetime(obs["time"].values).floor("s")
    n_points = len(times)
    rows = []
    for p in range(n_points):
        surface = None if surface_pressure_hpa is None else surface_pressure_hpa[p]
        centre_pa, thickness_pa = release_layers_pa(obs["pressure_levels"].values[p], surface)
        lat = float(obs["lat"].values[p])
        lon = float(obs["lon"].values[p])
        if box_degrees is None:
            dlon, dlat = footprint_box_degrees(lat, footprint_diameter_km)
        else:
            dlon, dlat = box_degrees[0][p], box_degrees[1][p]
        for level in range(1, n_levels + 1):
            rows.append(
                {
                    "source_name": f"SourceID_{p * n_levels + level}",
                    "point_number": p + 1,
                    "level_number": level,
                    "release_end": times[p],
                    "release_start": times[p] + pd.Timedelta(seconds=release_seconds),
                    "lon": lon,
                    "lat": lat,
                    "dlon": float(dlon),
                    "dlat": float(dlat),
                    "z_pa": centre_pa[level - 1],
                    "dz_pa": thickness_pa[level - 1],
                    "particles_per_second": particles_per_second(level, level_split),
                }
            )
    return pd.DataFrame(rows)
