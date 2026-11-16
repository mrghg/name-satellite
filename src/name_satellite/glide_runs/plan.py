"""Release plan for a GLIDE reference run: one NAME source per hourly release at one receptor.

Where this sits
---------------
    spec.py (receptors, blocks)  ->  **plan.py: release plan**  ->  prepare.py (NAME input, job)

The plan is a table with one row per release, in the columns ``name_input`` reads:

* ``source_name``: ``SourceID_<n>``, n = 1, 2, ... in time order. NAME's field request for it is
  ``Req_H1Z1_<n>``, so ``deliver.py`` matches footprints to releases by this name;
* ``release_time``: the release's name in the spec (§4): the START of its hour, UTC;
* ``release_start``, ``release_end``: NAME's source Start and Stop Time. In a backwards run a
  source "starts" at the later time, so ``release_start`` = release_time + 1 h and
  ``release_end`` = release_time (the same convention as the satellite plan);
* ``receptor``, ``lon``, ``lat`` (degrees);
* ``z_coord``, ``z``, ``dz``: vertical coordinate, centre and full depth. A point release (§3) is
  "m agl" with ``dz`` = 0;
* ``dlon``, ``dlat``: horizontal size, 0 for a point (§3);
* ``particles_per_second``: NAME's release rate (§4).
"""

import pandas as pd

from name_satellite.glide_runs import spec


def point_release_plan(receptor, release_times):
    """Release plan for a point receptor.

    Args:
        receptor: a ``spec.Receptor`` with a height (not the Lamont column).
        release_times: start times of the releases (UTC, on the hour); each release lasts
            ``spec.RELEASE_HOURS``.

    Returns:
        DataFrame, one row per release (columns in the module docstring).

    Raises:
        ValueError for the column receptor, an empty or unsorted list of times, repeated times or
        times not on the hour.
    """
    if receptor.is_column:
        raise ValueError(f"{receptor.code} is a column receptor; use the column plan (§5)")
    times = pd.DatetimeIndex(release_times)
    if len(times) == 0:
        raise ValueError("No release times given")
    if not times.is_monotonic_increasing or times.has_duplicates:
        raise ValueError("Release times must be increasing and distinct")
    if (times != times.floor("1h")).any():
        raise ValueError("Releases start on the hour (§4)")

    hour = pd.Timedelta(hours=spec.RELEASE_HOURS)
    n = len(times)
    return pd.DataFrame(
        {
            "source_name": [f"SourceID_{i}" for i in range(1, n + 1)],
            "release_time": times,
            "release_start": times + hour,  # NAME's Start Time: the later end of the hour
            "release_end": times,  # NAME's Stop Time: the earlier end
            "receptor": receptor.code,
            "lon": receptor.lon,
            "lat": receptor.lat,
            "z_coord": "m agl",
            "z": receptor.height_m_agl,
            "dz": 0.0,
            "dlon": 0.0,
            "dlat": 0.0,
            "particles_per_second": spec.POINT_PARTICLES_PER_SECOND,
        }
    )
