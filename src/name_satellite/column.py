"""Step 6c of the workflow: combine the per-level results into one column result per sounding.

Where this sits (see docs/workflow.md)
--------------------------------------
    6a. footprints, 6b. boundary histograms  ->  **6c. column combination**  ->  write output file

A satellite measures a weighted average of the mole fraction through the column. The weight of
retrieval level j is (pressure weight)_j × (averaging kernel)_j. So the column footprint (and
likewise the column boundary histograms) is

    column = Σ_j  w_j × (result of the NAME release that represents level j),   j < max_level

Levels from ``max_level`` upwards are left out; the inversion fills them from the a priori.

Legacy method for choosing "the NAME release that represents level j"
---------------------------------------------------------------------
Legacy (NAME_processing ``satellite_vertical_profile``; identical in acrg 75793919) picks, for
each retrieval level, the NAME release whose release-point pressure is NEAREST to the retrieval
level's pressure. The release-point pressure comes from the legacy met run, and equals the release
layer's centre unless the layer centre was below the model ground (docs/deviations.md X9). With
one release per retrieval level this is usually, but not always, the release for that same level;
see ``nearest_release_levels``.
"""

import numpy as np


def nearest_release_levels(release_pressure_pa, retrieval_pressure_hpa, max_level):
    """For each retrieval level below ``max_level``, the index of the nearest NAME release level.

    Args:
        release_pressure_pa: pressure (Pa) at each NAME release of one sounding (from the met run).
        retrieval_pressure_hpa: retrieval level pressures (hPa) of the sounding, surface first.
        max_level: number of retrieval levels to use (levels 0 .. max_level-1).

    Returns:
        Integer array of length ``max_level``: element j is the 0-based NAME release level used for
        retrieval level j. Ties go to the lower index (``numpy.argmin``), as in legacy.
    """
    release_pressure_pa = np.asarray(release_pressure_pa, dtype=float)
    retrieval_pa = np.asarray(retrieval_pressure_hpa, dtype=float)[:max_level] * 100.0
    return np.array([np.abs(release_pressure_pa - p).argmin() for p in retrieval_pa])


def column_weights(pressure_weights, averaging_kernel, max_level):
    """Weight of each retrieval level below ``max_level``: pressure_weight × averaging_kernel."""
    return np.asarray(pressure_weights)[:max_level] * np.asarray(averaging_kernel)[:max_level]


def combine_levels(per_level, chosen_levels, pressure_weights, averaging_kernel):
    """Weighted sum over retrieval levels of the chosen NAME releases' results.

    Args:
        per_level: array whose first dimension is the NAME release level (e.g. footprints with shape
            (n_release_levels, lat, lon)).
        chosen_levels: from ``nearest_release_levels`` (one entry per retrieval level used).
        pressure_weights, averaging_kernel: the sounding's retrieval values, surface first, as
            stored in the observation file.

    Returns:
        Σ_j per_level[chosen_levels[j]] × pressure_weight[j] × averaging_kernel[j], with the first
        dimension removed.
    """
    # Same order of operations as legacy -- (value × pressure weight) × averaging kernel, summed
    # from level 0 upwards -- so that floating-point rounding matches the legacy files exactly.
    total = per_level[chosen_levels[0]] * pressure_weights[0] * averaging_kernel[0]
    for j in range(1, len(chosen_levels)):
        total = total + per_level[chosen_levels[j]] * pressure_weights[j] * averaging_kernel[j]
    return total
