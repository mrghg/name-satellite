"""Step 6b, LEGACY method: where particles left the domain (boundary sensitivities).

Used by ``legacy_pipeline`` (validation only). New runs use the particle-fate accounting of
decision D13 in ``processing.py`` (``classify_fates``, ``release_fates``).

Where this sits (see docs/workflow.md)
--------------------------------------
    5. read NAME output  ->  **6b. boundary histograms**  ->  6c. combine into column footprints

The inversion estimates each observation's "background" by multiplying these histograms by mole
fractions from a global model at the domain edges. For each release, the histograms say what
fraction of its particles left through each edge (north, east, south, west), and where along the
edge and at what height.

This module implements only the LEGACY method, to reproduce the legacy files exactly.
Its known problems are described in docs/deviations.md (X6, X14). In short:

* a particle is assigned to an edge if its final position lies beyond the outermost half-cell of
  the OUTPUT grid on that side (north: lat > lat_max - dlat/2, etc.), so a particle near a corner
  can be counted on two edges;
* along each edge the particles are binned by grid cell and by height (20 bins, 0-20 km, m above
  ground) with ``numpy.histogram2d``, which silently drops particles outside the bins;
* each release's four histograms are divided by their combined count, so they sum to 1;
* if no particle reached any edge, the histograms of the previous level are copied (or, for the
  first level of a sounding, those of level 1 of the previous sounding).

Legacy source: NAME_processing ``particle_locations`` (today's version) and acrg 75793919 (the code
that made the 2019 production files). The only difference between them for this purpose is how the
far end of each edge is binned; choose it with ``variant``:

* ``"2019"``: the outermost bin ends at the grid edge (particles beyond it are dropped);
* ``"2025"``: the outermost bin is stretched to include all particles beyond the far end
  (NAME_processing commit 094d8e6, May 2025). Particles beyond the near end are still dropped.
"""

import numpy as np

# Legacy height bins for particles leaving through the lateral edges: centres 500, 1500, ... 19500 m
# (above ground level), i.e. 20 bins of 1 km from 0 to 20 km. Written exactly as in the legacy code.
LEGACY_HEIGHT_CENTRES_M = np.arange(0, 19001, 1000) + 1000 / 2.0


def _bin_edges(centres, values, variant):
    """Bin edges halfway between ``centres``, as in the legacy ``particle_location_edges``."""
    spacing = centres[1] - centres[0]
    edges = np.append(centres - spacing / 2.0, centres[-1] + spacing / 2.0)
    if variant == "2025" and values.size > 0:
        edges[-1] = max(np.max(values), edges[-1])
    elif variant not in ("2019", "2025"):
        raise ValueError(f"Unknown legacy variant {variant!r}; use '2019' or '2025'")
    return edges


def _edge_histogram(along_edge, height, along_edge_centres, height_centres, variant):
    """Count particles in (along-edge cell, height bin); returns shape (n_along, n_height)."""
    along_edges = _bin_edges(along_edge_centres, along_edge, variant)
    height_edges = _bin_edges(height_centres, height, variant)
    counts, _, _ = np.histogram2d(along_edge, height, bins=(along_edges, height_edges))
    return counts


def legacy_boundary_histograms(particles, lat, lon, n_sources, variant="2019"):
    """Legacy boundary-particle histograms for every release (NAME source) of one run.

    Args:
        particles: DataFrame from ``name_output.read_acrg_particle_file`` (columns ``lon``, ``lat``,
            ``height_m_agl``, ``source_number``).
        lat, lon: cell-centre coordinates of the output grid (degrees). For bit-exact legacy
            results build them with ``legacy_files.legacy_grid_coordinates``.
        n_sources: number of NAME sources in the run; sources are numbered 1..n_sources.
        variant: "2019" or "2025" (see module docstring).

    Returns:
        dict with
          * ``"north"``, ``"south"``: arrays (n_sources, n_lon, n_height);
          * ``"east"``, ``"west"``: arrays (n_sources, n_lat, n_height);
          * ``"height"``: the height-bin centres (m above ground);
          * ``"copied_from"``: array (n_sources,) giving, for each source whose particles never
            reached an edge, the source number whose histograms were copied (0 otherwise).
        Each source's four histograms sum to 1 (fractions of the particles counted on the edges).
    """
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    heights = LEGACY_HEIGHT_CENTRES_M
    dlat = lat[1] - lat[0]
    dlon = lon[1] - lon[0]

    # Which particles count for each edge (legacy test against the OUTPUT grid's outer half-cell).
    on_edge = {
        "north": particles["lat"] > lat.max() - dlat / 2.0,
        "east": particles["lon"] > lon.max() - dlon / 2.0,
        "south": particles["lat"] < lat.min() + dlat / 2.0,
        "west": particles["lon"] < lon.min() + dlon / 2.0,
    }
    # North/south edges are binned along longitude; east/west along latitude.
    along = {
        "north": ("lon", lon),
        "south": ("lon", lon),
        "east": ("lat", lat),
        "west": ("lat", lat),
    }

    result = {edge: np.zeros((n_sources, len(along[edge][1]), len(heights))) for edge in on_edge}
    copied_from = np.zeros(n_sources, dtype=int)

    for source_number in range(1, n_sources + 1):
        index = source_number - 1
        of_this_source = particles["source_number"] == source_number
        for edge, selected in on_edge.items():
            subset = particles[selected & of_this_source]
            column, centres = along[edge]
            result[edge][index] = _edge_histogram(
                subset[column].to_numpy(),
                subset["height_m_agl"].to_numpy(),
                centres,
                heights,
                variant,
            )

        total = sum(result[edge][index].sum() for edge in on_edge)
        if total > 0:
            for edge in on_edge:
                result[edge][index] = result[edge][index] / total
        elif source_number > 1:
            # No particle reached an edge: legacy copies another release's histograms
            # (see _legacy_copy_source for which one). docs/deviations.md X16.
            copy_index = _legacy_copy_source(source_number) - 1
            for edge in on_edge:
                result[edge][index] = result[edge][copy_index]
            copied_from[index] = copy_index + 1
        # A first source with no edge particles stays all zero (legacy behaviour).

    result["height"] = heights
    result["copied_from"] = copied_from
    return result


LEGACY_LEVELS_PER_SOUNDING = 17


def _legacy_copy_source(source_number, levels=LEGACY_LEVELS_PER_SOUNDING):
    """Source number whose histograms the legacy code copies when a source has no edge particles.

    Legacy (``particle_locations``): for level l > 1 copy level l-1 of the same sounding; for
    level 1 copy level 1 ("lev 0") of the previous sounding.
    """
    point_index, level_index = divmod(source_number - 1, levels)
    if level_index > 0:
        return source_number - 1
    return (point_index - 1) * levels + 1
