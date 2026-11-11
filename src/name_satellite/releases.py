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

Surface pressure: retrieval, model-scaled (D19) or legacy
---------------------------------------------------------
The retrieval's surface pressure p[0] and the met model's differ (smoothed model topography): on
1 Jan 2013 the model's was a median 10 hPa higher (model ground ~85 m lower), so layers taken
straight from the retrieval usually start above the model ground, and occasionally below it.

* Default for new runs (decision D19): pass ``model_surface_pressure_hpa`` (from the met-only NAME
  run, ``surface_pressure.py``). EVERY level is multiplied by model / retrieval surface pressure,
  so the layers follow the model's terrain: layer 0 starts at the model ground, and each layer
  keeps exactly its fraction of the column, so the retrieval's pressure weights and averaging
  kernel still apply unchanged. Every level moves by the same fraction of its pressure: median
  1%, up to 3.5% on 1 Jan 2013.
* Without it, the retrieval's own levels are used unchanged.
* Legacy (validation only): the 2013 legacy runs replaced p[0] alone by a surface pressure from a
  separate NAME run (``use_name_pressure=True``; docs/deviations.md X9). Pass
  ``surface_pressure_hpa`` to reproduce that.
"""

import numpy as np
import pandas as pd

HPA_TO_PA = 100.0


def scale_to_model_surface(pressure_levels_hpa, model_surface_pressure_hpa):
    """Retrieval levels scaled so that level 0 equals the model's surface pressure (D19).

    Args:
        pressure_levels_hpa: retrieval level pressures in hPa, surface first.
        model_surface_pressure_hpa: the met model's surface pressure at the sounding (hPa).

    Returns:
        New array: every level multiplied by model / retrieval surface pressure. Layer boundaries
        (``layer_boundaries_hpa``) scale by the same factor, so each layer keeps its fraction of
        the column.
    """
    p = np.array(pressure_levels_hpa, dtype=float)
    if not model_surface_pressure_hpa > 0:
        raise ValueError(f"Model surface pressure must be > 0; got {model_surface_pressure_hpa}")
    return p * (model_surface_pressure_hpa / p[0])


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


# Decision D21: NAME models the levels whose release layers lie wholly inside the dispersion domain
# (20 km agl top): 16 for the University of Leicester GOSAT-1 v10.0 and GOSAT-2 v2.0 retrievals
# (layer 16's top, about 65 hPa, is ~1 km below the domain top over Brazil). Levels above use the
# retrieval's prior in the inversion (Ganesan et al. 2017). The legacy runs released 17 (X10).
MODELLED_LEVELS = {"GOSAT-1 v10.0": 16, "GOSAT-2 v2.0": 16}
# Safety net for other retrievals: the modelled levels must carry at least this much of each
# sounding's column weight (pressure weight x averaging kernel). Jan 2013 GOSAT-1 v10: min 0.884,
# median 0.929; GOSAT-2 1 Jan 2020: min 0.887.
MIN_MODELLED_WEIGHT = 0.85


def check_modelled_levels(obs, scaled_levels_hpa, domain_top_hpa, n_levels,
                          min_weight=MIN_MODELLED_WEIGHT):  # fmt: skip
    """Check, per sounding, that the modelled levels are inside the domain and carry the column.

    D21: (1) the top of the highest modelled release layer (D19-scaled) must be below the domain
    top, i.e. at a higher pressure than the model pressure there; (2) the modelled levels must
    carry at least ``min_weight`` of the column weight, Σ pressure weight × averaging kernel.

    Args:
        obs: observation Dataset (``pressure_weights``, ``xch4_averaging_kernel``).
        scaled_levels_hpa: (sounding, level) release levels after D19 scaling (hPa).
        domain_top_hpa: model pressure at the domain top per sounding (hPa; met-only run).
        n_levels: number of modelled levels.
        min_weight: smallest acceptable modelled column weight.

    Raises:
        ValueError naming the soundings that fail.
    """
    ids = [str(i) for i in obs["sounding_id"].values]
    layer_top = np.array([layer_boundaries_hpa(levels)[n_levels] for levels in scaled_levels_hpa])
    outside = np.flatnonzero(layer_top <= np.asarray(domain_top_hpa))
    if outside.size:
        raise ValueError(
            f"Release layer {n_levels} reaches above the domain top for soundings "
            f"{[ids[i] for i in outside]}: layer top {layer_top[outside]} hPa, domain top "
            f"{np.asarray(domain_top_hpa)[outside]} hPa (raise the domain top or model fewer "
            "levels)"
        )
    weight = obs["pressure_weights"].values * obs["xch4_averaging_kernel"].values
    modelled = weight[:, :n_levels].sum(axis=1)
    low = np.flatnonzero(modelled < min_weight)
    if low.size:
        raise ValueError(
            f"The {n_levels} modelled levels carry less than {min_weight} of the column weight for "
            f"soundings {[ids[i] for i in low]}: {modelled[low]} (check the retrieval's levels)"
        )


# Decision D22: particle release rate (particles/s) of each modelled level, per product. Levels
# 1-13 release 30/s, level 14 10/s and levels 15-16 3/s (legacy: 30/s up to level 12, 3/s above,
# ``particles_per_second``). On the 1 Jan 2013 seed pair, level 13 (around 170 hPa) at 3/s gave
# ~40% of the footprint noise and ~36% of the boundary noise for ~4% of the footprint; level 14
# mainly adds boundary noise; levels 15-16 add little to either. This allocation cuts the noise
# variance by ~37% (footprint) and ~40% (boundary) for ~10% more particle-hours.
PARTICLE_RATES = {
    "GOSAT-1 v10.0": (30,) * 13 + (10, 3, 3),
    "GOSAT-2 v2.0": (30,) * 13 + (10, 3, 3),
}


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
    model_surface_pressure_hpa=None,
    particle_rates=None,
):
    """One row per NAME release (source) for one day of soundings.

    Args:
        obs: observation Dataset, one entry per sounding, in time order: our format
            (``observations.to_obs``, dim ``sounding_id``) or the legacy one (dim ``time``). Uses
            ``time``, ``lon``, ``lat`` and ``pressure_levels`` (hPa).
        n_levels: retrieval levels released per sounding, from the surface (legacy 17).
        release_seconds: length of each release (legacy 60 s).
        level_split, footprint_diameter_km: see ``particles_per_second``, ``footprint_box_degrees``.
        particle_rates: optional particles/s for each level 1..n_levels (e.g.
            ``PARTICLE_RATES[product]``, D22); default the legacy ``level_split`` rule.
        surface_pressure_hpa: optional per-sounding replacement of level 0 (legacy X9;
            validation runs only).
        model_surface_pressure_hpa: optional per-sounding model surface pressure (hPa); every
            level is scaled to it (``scale_to_model_surface``, D19). The default for new runs.
        box_degrees: optional (dlon, dlat) arrays per sounding, overriding the footprint box (e.g.
            the widened boxes of legacy binned soundings).

    Returns:
        DataFrame with columns
          * ``source_name``: "SourceID_<n>", n = 1, 2, ... in sounding order, levels within each;
          * ``point_number``, ``level_number``: 1-based sounding and level numbers;
          * ``sounding_id``: the sounding's ID (D18), only if ``obs`` has one;
          * ``release_end``: the sounding time truncated to the whole second (legacy, UTC). In
            NAME's backwards run this is the source's "Stop Time";
          * ``release_start``: ``release_end`` + ``release_seconds`` (NAME's "Start Time");
          * ``lon``, ``lat``, ``dlon``, ``dlat``: box centre and full widths (degrees);
          * ``z_pa``, ``dz_pa``: layer centre and thickness (Pa);
          * ``particles_per_second``.
    """
    times = pd.to_datetime(obs["time"].values).floor("s")
    n_points = len(times)
    sounding_ids = obs["sounding_id"].values if "sounding_id" in obs.variables else None
    if surface_pressure_hpa is not None and model_surface_pressure_hpa is not None:
        raise ValueError("Give surface_pressure_hpa (legacy) or model_surface_pressure_hpa")
    if particle_rates is None:
        rates = [particles_per_second(level, level_split) for level in range(1, n_levels + 1)]
    else:
        rates = [int(rate) for rate in particle_rates]
        if len(rates) != n_levels or min(rates) < 1:
            raise ValueError(f"Need {n_levels} positive particle rates; got {particle_rates}")
    rows = []
    for p in range(n_points):
        levels = obs["pressure_levels"].values[p]
        if model_surface_pressure_hpa is not None:
            levels = scale_to_model_surface(levels, model_surface_pressure_hpa[p])
        surface = None if surface_pressure_hpa is None else surface_pressure_hpa[p]
        centre_pa, thickness_pa = release_layers_pa(levels, surface)
        lat = float(obs["lat"].values[p])
        lon = float(obs["lon"].values[p])
        if box_degrees is None:
            dlon, dlat = footprint_box_degrees(lat, footprint_diameter_km)
        else:
            dlon, dlat = box_degrees[0][p], box_degrees[1][p]
        for level in range(1, n_levels + 1):
            row = {
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
                "particles_per_second": rates[level - 1],
            }
            if sounding_ids is not None:
                row["sounding_id"] = str(sounding_ids[p])
            rows.append(row)
    return pd.DataFrame(rows)
