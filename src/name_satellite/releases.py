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
