"""Measures for comparing two sets of processed footprints (validation of new NAME runs).

Where this sits (see docs/workflow.md)
--------------------------------------
VALIDATION ONLY, not part of the production workflow. Used after step 6 (by the ``dev/scripts``
stage scripts) to compare footprints from two runs, e.g.
NAME 8.7.1 against the legacy production files, or two runs that differ only in their random seed
(which measures pure particle noise).

NAME is a Monte Carlo model, so two runs can only be compared statistically. The question is
always: is run X - reference no larger than the difference between two runs that differ only by
chance? Hence every measure here is applied to both pairs and the results set side by side.

All functions take footprints as arrays with dims (sounding, lat, lon), in any consistent units.
"""

import numpy as np


def total_ratio(f, g):
    """Per-sounding ratio of footprint totals, Σf / Σg (array, one value per sounding)."""
    return f.sum(axis=(1, 2)) / g.sum(axis=(1, 2))


def near_field_ratio(f, g, lat, lon, sounding_lat, sounding_lon, half_width_deg=2.0):
    """Per-sounding ratio of footprint totals within ±half_width_deg of each sounding.

    Args:
        f, g: footprints (sounding, lat, lon).
        lat, lon: grid coordinates (degrees).
        sounding_lat, sounding_lon: position of each sounding (degrees).
    """
    ratios = []
    for p in range(f.shape[0]):
        near = (np.abs(lat[:, None] - sounding_lat[p]) < half_width_deg) & (
            np.abs(lon[None, :] - sounding_lon[p]) < half_width_deg
        )
        ratios.append(f[p][near].sum() / g[p][near].sum())
    return np.array(ratios)


def coarsen(f, factor):
    """Sum blocks of factor × factor grid cells (cells left over at the edges are dropped)."""
    n, ny, nx = f.shape
    ny2, nx2 = ny // factor * factor, nx // factor * factor
    return f[:, :ny2, :nx2].reshape(n, ny2 // factor, factor, nx2 // factor, factor).sum((2, 4))


def shape_difference(f, g, factor=1):
    """Per-sounding L1 distance between the footprint SHAPES of f and g.

    Each footprint is first divided by its own total (so differences in totals don't count), and
    optionally coarsened by ``factor``. Then Σ|f - g| over the grid: 0 = identical, 2 = no
    overlap. Particle noise falls as the grid is coarsened; a systematic difference doesn't.
    """
    f, g = coarsen(f, factor), coarsen(g, factor)
    f = f / f.sum(axis=(1, 2), keepdims=True)
    g = g / g.sum(axis=(1, 2), keepdims=True)
    return np.abs(f - g).sum(axis=(1, 2))


def histogram_difference(a, b, names=("n", "e", "s", "w")):
    """Per-sounding Σ|a - b| over the four edge histograms, relative to Σ b.

    Args:
        a, b: dicts {edge: array (sounding, position along edge, height)}.
    """
    diff = sum(np.abs(a[n] - b[n]).sum(axis=(1, 2)) for n in names)
    total = sum(b[n].sum(axis=(1, 2)) for n in names)
    return diff / total
