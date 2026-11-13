"""Met at the release points, from the main NAME run (replaces the legacy separate met run).

Where this sits (see docs/workflow.md)
--------------------------------------
    3. write NAME input  (**release-met output requests added here**)  ->  4. run NAME
    5. read NAME output  (**release-met read here**)  ->  6. processing

What is needed, and why
-----------------------
For each release point (in new runs: each release group, D20; the rows are still called
"soundings" here), the processing needs:

* the model's surface pressure at the sounding: a release layer centred below it is partly
  underground (deviation X9), and NAME then reports a pressure for the release other than the
  one requested, which the legacy column combination used to pick levels. (New runs get it
  beforehand too, from the met-only run in ``surface_pressure.py``, which writes and reads its
  requests with this module; the two should agree, D19);
* temperature, pressure, wind speed and direction at the lowest release, and the boundary-layer
  depth: stored in the footprint file (legacy kept those of the lowest release, level 1).

The legacy method was a second NAME run with one horizontal, vertical and time grid per release
(Satellite_launch.scr). NAME 8.7.1 allows only 70/20/200 such grids, so it fails for any busy
day. Here instead, in the main run:

* the soundings are split into groups of at most MAX_PER_GROUP (NAME's longest input array);
  each group has ONE set of locations (the sounding positions) and ONE vertical grid (an array
  of the level-1 release pressures, Pa), shared by its soundings;
* each sounding has its own time grid (one time) and output group (one file). NAME 8.7.1 can't
  output fields on an irregular time array (a T-Array grid aborts with "UNEXPECTED FATAL ERROR in
  MultiplyShortTime"), so one time grid per sounding it is;
* NAME outputs every combination of location and pressure at each sounding's time; we keep only
  the sounding's own (location, pressure). Cost: n² values per quantity and sounding for a group
  of n soundings, n³ in all, at most 86³ ≈ 640,000: small.

Limits per run (NAME 8.7.1): 200 time grids (one is the footprints' own), so MAX_SOUNDINGS =
199; 5 sets of locations and 12 input arrays (one per group),
which allow more. 199 soundings is about a week of GOSAT-BRAZIL, matching the 3,000-source limit
of batching.py. The output groups (one per sounding) need ``RunSettings.max_field_output_groups``.

Times: the met is sampled at the sounding time truncated to the whole second (``release_end``),
as legacy (TGrid t0 = sounding time).

Output files
------------
One NAME point-output file per release point, ``<output folder>/ReleaseMet_<S0001...>_C1.txt``
(S0001 = first row of ``release_met_soundings``, in release-plan order). In each, the header rows
give per column the request (``<S...>_T``, ``_P``, ``_WS``, ``_WD``, ``_BL``, ``_PS``), quantity,
units and "Z = ..." level; then one data row per location of the output group. Only the row of
the file's own location, and for the level-1 quantities the column at its own release pressure,
are used (``read_release_met``). ``processing.read_day`` reads the main run's files; the met-only
run (``surface_pressure``) writes the same kind of file with ``_PS`` and ``_PT`` only.
"""

from pathlib import Path

import numpy as np
import pandas as pd

MAX_PER_GROUP = 86  # NAME 8.7.1: longest input array
MAX_SOUNDINGS = 199  # NAME 8.7.1: 200 time grids, less the footprints' own
NAME_TIME_FORMAT = "%d/%m/%Y %H:%M:%S"

# (request suffix, NAME quantity, vertical grid: "release" = level-1 pressure, "ground" = 0 m agl,
# "top" = the dispersion domain's top, m agl: only with ``top_m_agl``)
QUANTITIES = [
    ("T", "Temperature (K)", "release"),
    ("P", "Pressure (Pa)", "release"),
    ("WS", "Wind Speed", "release"),
    ("WD", "Wind Direction (degrees)", "release"),
    ("BL", "Boundary Layer Depth", "ground"),
    ("PS", "Pressure (Pa)", "ground"),
]


def release_met_soundings(plan):
    """One row per sounding: where and when to sample met.

    Args:
        plan: release plan (``releases.build_release_plan``, or a combined plan from
            ``batching.combine_release_plans``, which has a ``day`` column).

    Returns:
        DataFrame with ``sounding`` (a unique label, e.g. "S0001", in plan order), ``day`` (if in
        the plan), ``point_number``, ``release_group`` or ``sounding_id`` (whichever the plan has),
        ``lon``, ``lat``
        (degrees), ``time`` (sounding time, truncated to the second) and ``z_pa`` (centre of
        the level-1 release layer, Pa), plus ``group`` (1-based) assigning soundings to output
        groups of at most MAX_PER_GROUP.

    Raises:
        ValueError if there are more than MAX_SOUNDINGS soundings.
    """
    lowest = plan[plan["level_number"] == 1].reset_index(drop=True)
    columns = ["day", "point_number"] if "day" in lowest else ["point_number"]
    # Carry the plan's ID: release_group (new runs) or sounding_id (per-sounding plans).
    columns += [name for name in ("release_group", "sounding_id") if name in lowest]
    soundings = lowest[columns].copy()
    soundings.insert(0, "sounding", [f"S{i:04d}" for i in range(1, len(lowest) + 1)])
    soundings["lon"] = lowest["lon"]
    soundings["lat"] = lowest["lat"]
    soundings["time"] = lowest["release_end"]
    soundings["z_pa"] = lowest["z_pa"]
    soundings["group"] = soundings.index // MAX_PER_GROUP + 1
    if len(soundings) > MAX_SOUNDINGS:
        raise ValueError(
            f"{len(soundings)} soundings in one run; release met allows at most {MAX_SOUNDINGS}"
        )
    return soundings


def release_met_blocks(soundings, dx, dy, quantities=QUANTITIES, top_m_agl=None):
    """NAME input blocks requesting met at each sounding (text to add to the input file).

    Args:
        soundings: from ``release_met_soundings``.
        dx, dy: nominal grid spacing for the location grids (degrees). NAME requires one for
            grids built from locations; it doesn't affect point met.
        quantities: which met to request (default QUANTITIES; the met-only surface-pressure run
            asks for surface pressure only, ``surface_pressure.py``).
        top_m_agl: height (m agl) of the "top" vertical grid: the dispersion domain's top, for
            quantities with level "top" (D21).

    Output: one file per sounding, output group ``ReleaseMet_<sounding>``.
    """
    at_release = any(level == "release" for _, _, level in quantities)
    at_top = any(level == "top" for _, _, level in quantities)
    if at_top and top_m_agl is None:
        raise ValueError("A quantity at the domain top needs top_m_agl")
    lines = [
        "! Met at the release points (name_satellite.release_met)",
        "Vertical Grids:",
        "Name, Z-Coord, nZ, dZ, Z0",
        "ZGrid_Ground, m agl, 1, 0.0, 0.0",
    ]
    if at_top:
        lines.append(f"ZGrid_Top, m agl, 1, 0.0, {float(top_m_agl)!r}")
    lines.append("")
    for group, members in soundings.groupby("group"):
        g = f"RM{group}"
        lines += [f"Locations: ReleasePoints_{group}", "Name, H-Coord, X, Y"]
        lines += [f"{s.sounding}, Lat-Long, {s.lon!r}, {s.lat!r}" for s in members.itertuples()]
        lines += [
            "",
            "Horizontal Grids:",
            "Name, dX, dY, Set of Locations",
            f"HGrid_{g}, {dx}, {dy}, ReleasePoints_{group}",
            "",
        ]
        if at_release:
            # A NAME grid is a set of values: sort the pressures (surface first), drop repeats.
            pressures = sorted(set(members["z_pa"]), reverse=True)
            lines += [
                f"Array: ZArray_{g}",
                "Array Values",
                *[repr(p) for p in pressures],
                "",
                "Vertical Grids:",
                "Name, Z-Coord, Z-Array",
                f"ZGrid_{g}, Pa, ZArray_{g}",
                "",
            ]
        lines += [
            "Temporal Grids:",
            "Name, nT, dT, T0",
            *[
                f"TGrid_{s.sounding}, 1, 00:01, {s.time.strftime(NAME_TIME_FORMAT)}"
                for s in members.itertuples()
            ],
            "",
            "Output Requirements - Fields:",
            "Name, Quantity, Species, Source, H-Coord, H-Grid, Z-Grid, T-Grid, BL Average, "
            "T Av Or Int, Av Time, # Av Times, Sync?, Output Route, Across, Separate File, "
            "Output Format, Output Group",
        ]
        for s in members.itertuples():
            for suffix, quantity, level in quantities:
                z_grids = {"release": f"ZGrid_{g}", "ground": "ZGrid_Ground", "top": "ZGrid_Top"}
                z_grid = z_grids[level]
                # Wind direction needs a horizontal coordinate system to be measured in.
                h_coord = "Lat-Long" if "Direction" in quantity else ""
                lines.append(
                    f"{s.sounding}_{suffix}, {quantity}, , , {h_coord}, HGrid_{g}, {z_grid}, "
                    f"TGrid_{s.sounding}, No, No, , , Yes, D, Z, , IA, ReleaseMet_{s.sounding}"
                )
        lines.append("")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------------------------
# Reading the output
# ---------------------------------------------------------------------------------------------

# Columns written by ``release_met_blocks``, and the names we give them.
OUTPUT_NAMES = {
    "T": "temperature_k",
    "P": "pressure_pa",
    "WS": "wind_speed_m_s",
    "WD": "wind_direction_deg",
    "BL": "boundary_layer_depth_m",
    "PS": "surface_pressure_pa",
    "PT": "domain_top_pressure_pa",
}


def read_release_met_file(path):
    """Read one release-met output file (NAME fields format on a set of locations).

    Returns:
        (columns, values): ``columns`` is a DataFrame with one row per field column (``request``,
        ``quantity``, ``units``, ``z`` as printed, e.g. "Z = 99886.84 Pa"); ``values`` is a
        DataFrame indexed by location name with one column per field column (same order).
    """
    with open(path) as f:
        lines = f.read().splitlines()
    fields_line = lines.index(" Fields:")
    heading = next(i for i, line in enumerate(lines) if line.lstrip().startswith("T Index"))
    rows = [[c.strip() for c in line.split(",")] for line in lines[fields_line + 1 : heading]]
    first = 4  # T Index, Location Index, T, Location
    z_rows = [r for r in rows if any(c.startswith("Z = ") for c in r)]
    if len(z_rows) != 1:
        raise ValueError(f"{path}: expected one 'Z = ' header row, found {len(z_rows)}")
    columns = pd.DataFrame(
        {
            "request": rows[1][first:],
            "quantity": rows[2][first:],
            "units": rows[4][first:],
            "z": z_rows[0][first:],
        }
    )
    columns = columns[columns["request"] != ""].reset_index(drop=True)  # trailing comma

    data = [[c.strip() for c in line.split(",")] for line in lines[heading + 1 :] if line.strip()]
    names = [row[3] for row in data]
    values = pd.DataFrame(
        [[float(c) for c in row[first : first + len(columns)]] for row in data], index=names
    )
    if values.index.duplicated().any():
        raise ValueError(f"{path}: a location appears more than once")
    return columns, values


def _printed_pressure(z_label):
    """Pressure in Pa from a header label like "Z = 99886.84 Pa"."""
    number, units = z_label.removeprefix("Z = ").split()
    if units != "Pa":
        raise ValueError(f"Expected a pressure level, got {z_label!r}")
    return float(number)


def read_release_met(output_folder, soundings, quantities=QUANTITIES):
    """Met at each sounding's own position, level-1 pressure and time.

    Args:
        output_folder: the NAME run's output folder (files ``ReleaseMet_<sounding>_C1.txt``).
        soundings: the table from ``release_met_soundings`` used to write the run.
        quantities: the quantities the run requested (as for ``release_met_blocks``).

    Returns:
        ``soundings`` with the columns in OUTPUT_NAMES added: temperature (K), pressure (Pa) at
        the level-1 release, wind speed (m/s) and direction (degrees), boundary-layer depth (m)
        and surface pressure (Pa), or those of ``quantities``. NAME prints 7 significant
        figures.

    Raises:
        FileNotFoundError if a file is missing; ValueError if a request or level is missing, or
        the level is ambiguous.
    """
    output_folder = Path(output_folder)
    results = []
    for s in soundings.itertuples():
        columns, values = read_release_met_file(output_folder / f"ReleaseMet_{s.sounding}_C1.txt")
        row = {}
        for suffix, _, _ in quantities:
            name = OUTPUT_NAMES[suffix]
            of_request = columns.index[columns["request"] == f"{s.sounding}_{suffix}"]
            if len(of_request) == 0:
                raise ValueError(f"{s.sounding}: no {s.sounding}_{suffix} columns")
            if len(of_request) == 1:  # ground-level quantities: a single column
                column = of_request[0]
            else:  # one column per pressure in the group: take this sounding's own
                printed = columns.loc[of_request, "z"].map(_printed_pressure)
                distance = (printed - s.z_pa).abs()
                column = distance.idxmin()
                # 7 significant figures: within 0.05 Pa; no other level that close.
                if distance.min() > 0.06 or (distance < 0.06).sum() > 1:
                    raise ValueError(f"{s.sounding}: can't identify level {s.z_pa} Pa")
            row[name] = values.loc[s.sounding, column]
        results.append(row)
    return pd.concat([soundings.reset_index(drop=True), pd.DataFrame(results)], axis=1)


def release_met_per_source(plan, release_met):
    """Release-point met for every source, in the layout of the legacy met files.

    This is what ``legacy_pipeline.process_legacy_inputs`` takes as ``met`` (as from
    ``legacy_files.read_legacy_met_folder``), so NAME 8.7.1 runs processed the legacy way need no
    separate met run. (``processing`` uses ``read_release_met`` directly instead.)

    Pressure: NAME reports the requested pressure for a release above the ground, and the model's
    surface pressure for one centred below it (X9). So for every release,
    ``pressure_pa = min(z_pa, surface pressure)``. Checked against the legacy met run on 1-2 Jan
    2013 (dev/scripts/legacy_validation/compare_release_met.py): equal to 0.02 Pa above ground;
    for the three underground releases of 2 Jan, NAME 8.7.1's pressure equals its surface pressure
    exactly and the legacy value within 2 Pa.

    Temperature, wind and boundary-layer depth are sampled at level 1 only (the legacy processing
    keeps only level 1's); other levels get NaN.

    Args:
        plan: one day's release plan (``releases.build_release_plan``), sources in order.
        release_met: ``read_release_met`` output for the same day (one row per sounding).

    Returns:
        DataFrame, one row per source: ``point_number``, ``level_number``, ``time``, ``lon``,
        ``lat``, ``z_pa``, ``temperature_c`` (°C), ``boundary_layer_depth_m``, ``pressure_pa``,
        ``wind_speed_m_s``, ``wind_direction_deg``.

    Raises:
        ValueError if the release met and the plan have different point numbers.
    """
    met = release_met.set_index("point_number")
    if sorted(met.index) != sorted(plan["point_number"].unique()):
        raise ValueError("Release met and release plan have different soundings")
    rows = []
    for release in plan.itertuples():
        sounding = met.loc[release.point_number]
        level_1 = release.level_number == 1
        rows.append(
            {
                "point_number": release.point_number,
                "level_number": release.level_number,
                "time": release.release_end,
                "lon": release.lon,
                "lat": release.lat,
                "z_pa": release.z_pa,
                "temperature_c": sounding.temperature_k - 273.15 if level_1 else np.nan,
                "boundary_layer_depth_m": sounding.boundary_layer_depth_m,
                "pressure_pa": min(release.z_pa, sounding.surface_pressure_pa),
                "wind_speed_m_s": sounding.wind_speed_m_s if level_1 else np.nan,
                "wind_direction_deg": sounding.wind_direction_deg if level_1 else np.nan,
            }
        )
    return pd.DataFrame(rows)
