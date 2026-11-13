"""Step 6a of the workflow: turn NAME fields into footprints, one per release.

Where this sits (see docs/workflow.md)
--------------------------------------
    5. read NAME output  ->  **6a. fields to footprints**  ->  6c. combine into column footprints

No files: takes the fields read by ``name_output`` (ppm s) and returns arrays, used by
``processing`` (new runs) and ``legacy_pipeline`` (legacy reproduction).

NAME releases particles at 1 g/s for a short time (60 s, legacy and new runs) and reports,
for each output grid cell, the mixing ratio in the near-surface layer (0-40 m) *integrated over
time*, in ppm s. Because the run goes backwards in time, that integral is the sensitivity of the
mixing ratio at the release point to a surface emission in that grid cell.

The footprint (the "source-receptor relationship") in (mol/mol) / (mol/m2/s) is

    footprint = field [ppm s] × 1e-6 [mol/mol per ppm] × cell area [m2] / release duration [s]

which assumes a release of 1 g/s and a molar mass of 1 g/mol (the tracer is defined that way in the
NAME runs). Legacy source: NAME_processing ``footprint_array.convert_units`` ('ppm s' branch), which
is identical in the code that made the 2019 production files (acrg 75793919).

Only mixing-ratio output is supported (decision D6): the legacy conversion of air concentration
(g s/m3) needs an air density, and the legacy choices for it were wrong (docs/deviations.md X15).
"""

import numpy as np

EARTH_RADIUS_M = 6367500.0  # value used by the legacy code (acrg_grid / NAME_processing grid.py)


def grid_cell_area_m2(lat, lon):
    """Area (m²) of each cell of a regular latitude-longitude grid.

    Identical to the legacy ``grid.areagrid``: each cell is the spherical band between lat ± dlat/2,
    times the longitude width, on a sphere of radius 6,367,500 m. The grid spacing is the mean
    spacing of the input coordinates.

    Args:
        lat: cell-centre latitudes (degrees), 1-D.
        lon: cell-centre longitudes (degrees), 1-D.

    Returns:
        2-D array of areas, shape (len(lat), len(lon)).
    """
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    dlon = abs(np.mean(lon[1:] - lon[:-1])) * np.pi / 180.0
    dlat = abs(np.mean(lat[1:] - lat[:-1])) * np.pi / 180.0
    colatitude = np.pi * (90.0 - lat) / 180.0

    area = np.zeros((len(lat), len(lon)))
    for i in range(len(lat)):
        if colatitude[i] == 0.0 or np.isclose(colatitude[i], np.pi):
            # A cell centred on a pole is a cap of half-width dlat/2.
            area[i, :] = EARTH_RADIUS_M**2 * abs(np.cos(dlat / 2.0) - np.cos(0.0)) * dlon
        else:
            north_edge = colatitude[i] - dlat / 2.0
            south_edge = colatitude[i] + dlat / 2.0
            area[i, :] = EARTH_RADIUS_M**2 * (np.cos(north_edge) - np.cos(south_edge)) * dlon
    return area


def footprint_from_mixing_ratio(field_ppm_s, area_m2, release_duration_s):
    """Convert time-integrated NAME mixing ratio into a footprint.

    Args:
        field_ppm_s: NAME field(s) in ppm s; the last two dimensions are (lat, lon).
        area_m2: grid-cell areas (m²), shape (lat, lon), from ``grid_cell_area_m2``.
        release_duration_s: length of the release (s), 60 for the legacy satellite runs.

    Returns:
        Footprint(s) in (mol/mol) / (mol/m2/s), same shape as ``field_ppm_s``.
    """
    # Same order of operations as the legacy code, so results agree to the last bit.
    return field_ppm_s * area_m2 * 1e-6 / release_duration_s
