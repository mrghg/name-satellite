"""Step 3 of the workflow: write the NAME input file for one day of satellite back-runs.

Where this sits (see docs/workflow.md)
--------------------------------------
    2. plan releases (releases.build_release_plan)  ->  **3. write NAME input**  ->  4. run NAME

The input file is plain text made of "blocks" (``Main Options:``, ``Sources:``, ...). NAME reads
the columns of each block by their heading, so the order of columns doesn't matter. This module
builds each block with one small function, then ``write_name_input`` joins them.

Syntax follows NAME 8.7.1 (the NAME team's 8.7.1 demo templates on the NAME GWS); the scientific
settings follow the legacy satellite runs unless a docs/deviations.md entry says otherwise:

* one source per release (cuboid, Lat-Long × Pa), INERT-TRACER at 1 g/s;
* one output request per source: mixing ratio in the 0-40 m layer, integrated over the whole run,
  on the output grid (one column per source in the fields file);
* particle fate output (where/when each particle stopped), replacing the ACRG-modified NAME's
  particle file;
* the UM global met declarations from the NAME team's 8.7.1 demo files.

What's different from the 2019 legacy input (deliberately):
* The time-integration window covers every release of the day (X18). Legacy anchored a 721 h
  window at the FIRST sounding, which seems to cut off the start of later soundings' trajectories.
  ``integration="legacy_2019"`` reproduces the old window.
"""

import math
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from name_satellite.met import DEMO_UTILITIES, NAME_DIR, NAME_GWS, MetVersion
from name_satellite.release_met import release_met_blocks, release_met_soundings

NAME_TIME_FORMAT = "%d/%m/%Y %H:%M:%S"


@dataclass(frozen=True)
class OutputGrid:
    """Horizontal output grid. ``x_min``/``y_min`` are the centre of the first cell (degrees)."""

    nx: int
    ny: int
    dx: float
    dy: float
    x_min: float
    y_min: float


@dataclass(frozen=True)
class DispersionDomain:
    """Region in which particles are followed; they stop when they leave it (degrees, m agl)."""

    x_min: float
    x_max: float
    y_min: float
    y_max: float
    z_max_m_agl: float = 19000.0


@dataclass(frozen=True)
class RunSettings:
    """Everything about a day's NAME run apart from the releases themselves.

    Defaults are the legacy satellite settings (2013 GOSAT-BRAZIL run) unless noted.
    """

    run_name: str
    output_folder: str
    grid: OutputGrid
    domain: DispersionDomain
    met_version: MetVersion
    met_folder: str  # where met files are restored to (scratch)
    max_age_hours: int = 720  # 30 days
    sync_minutes: int = 6
    max_particles: int = 7_000_000
    # NAME counts internal fields too (e.g. air density for mixing ratio, time integration), so
    # this must be well above the number of requests. Legacy value.
    max_field_requirements: int = 20_000
    # One output group for the footprints plus one per sounding for the release-point met
    # (release_met.py); NAME's default is 100.
    max_field_output_groups: int = 1000
    threads: int = 1
    random_seed: str = "Fixed"  # legacy used "Random (Parallel)"; "Fixed" makes runs repeatable
    min_bl_depth_m: float = 40.0  # legacy; the NAME demo declarations use 100 m
    surface_layer_depth_m: float = 40.0  # footprint layer 0-40 m (legacy ZGrid1)
    fate_z_coord: str = "m agl"  # legacy; decision D14 will switch to "Pa"
    # Add NAME's met columns (incl. Pressure, Pa) to each fate record: with fate_z_coord "m agl"
    # this gives both the height (to find top exits) and the pressure (D14 bins).
    fate_met: bool = False
    integration: str = "all_releases"  # or "legacy_2019" (X18)
    # Request met at the release points in this run (release_met.py). Legacy used a separate run.
    release_met: bool = True
    use_timers: bool = True  # NAME prints a summary of time spent in each part of the run
    topography_folder: str = str(NAME_GWS / "data" / "UMTopogData")
    met_restore_script: str = str(DEMO_UTILITIES / "MetRestore_JASMIN.sh")
    name_dir: str = str(NAME_DIR)


# ---------------------------------------------------------------------------------------------
# Timing
# ---------------------------------------------------------------------------------------------


def run_timing(plan, settings):
    """Key times of a day's run.

    In a backwards run, time runs from the latest release back to the earliest time any particle
    can reach. NAME integrates each field over one window ending ("end" in run time, i.e.
    chronologically earliest) at ``window_start``.

    Returns:
        dict with ``window_start`` (chronologically earliest time of the window, = earliest release
        end − max age), ``window_hours`` (length of the integration window, whole hours), and
        ``max_travel_hours`` (particles are stopped after this travel time).
    """
    earliest_end = plan["release_end"].min()
    latest_start = plan["release_start"].max()
    window_start = earliest_end - pd.Timedelta(hours=settings.max_age_hours)
    if settings.integration == "all_releases":
        # Cover from window_start up to the latest release start, rounded up to whole hours.
        span = (latest_start - window_start) / pd.Timedelta(hours=1)
        window_hours = math.ceil(span)
        max_travel_hours = settings.max_age_hours
    elif settings.integration == "legacy_2019":
        # Legacy (2019 runs): max age + 1 h from the first sounding, whatever the day's span (X18).
        window_hours = settings.max_age_hours + 1
        max_travel_hours = settings.max_age_hours + 1
    else:
        raise ValueError(f"Unknown integration option {settings.integration!r}")
    return {
        "window_start": window_start,
        "window_hours": window_hours,
        "max_travel_hours": max_travel_hours,
    }


# ---------------------------------------------------------------------------------------------
# Blocks
# ---------------------------------------------------------------------------------------------


def main_options_block(settings, n_sources):
    return (
        "Main Options:\n"
        "Absolute or Relative Time?, Fixed Met?, Flat Earth?, Run Name, Random Seed, "
        "Max # Sources, Max # Field Reqs, Max # Field Output Groups, Backwards?\n"
        f"Absolute, No, No, {settings.run_name}, {settings.random_seed}, "
        f"{n_sources}, {settings.max_field_requirements}, {settings.max_field_output_groups}, Yes\n"
    )


def restart_and_cases_block():
    # Required by NAME even when unused: no restart files, a single case (no ensembles). Copied
    # from the NAME team's 8.7.1 demo template. NAME needs the trailing commas: a line whose last
    # item is blank must end with a comma.
    return (
        "Restart File Options:\n"
        "# Cases Between Writes, Time Between Writes, Delete Old Files?, Write On Suspend?,\n"
        ", , , ,\n"
        "\n"
        "Multiple Case Options:\n"
        "Dispersion Options Ensemble Size, Met Ensemble Size,\n"
        "1, 1,\n"
    )


def openmp_block(settings):
    # "Threads" applies to every parallel loop (particles, particle update, output processing).
    # "Output Group Threads" parallelises over output groups; we have one, so 1 (as the NAME
    # demo). The job must provide the same number of CPUs and OMP_NUM_THREADS (lotus.py).
    use = "Yes" if settings.threads > 1 else "No"
    return (
        "OpenMP Options:\n"
        "Use OpenMP?, Threads, Output Group Threads\n"
        f"{use}, {settings.threads}, 1\n"
    )


def timer_block(settings):
    use = "Yes" if settings.use_timers else "No"
    return f"Timer Options:\nUse Timers?, Summary Only?\n{use}, Yes\n"


def files_and_coordinates_block(settings):
    definition = (
        Path(settings.name_dir) / "Resources" / "Defns" / settings.met_version.definition_file
    )
    return (
        "Output Options:\n"
        "Folder\n"
        f"{settings.output_folder}\n"
        "\n"
        "Input Files:\n"
        "File names\n"
        f"{definition}\n"
        "\n"
        "Horizontal Coordinate Systems:\n"
        "Name\n"
        "Lat-Long\n"
        "\n"
        "Vertical Coordinate Systems:\n"
        "Name\n"
        "m agl\n"
        "m asl\n"
        "Pa\n"
    )


def grids_block(settings, timing):
    g = settings.grid
    layer = settings.surface_layer_depth_m
    t0 = timing["window_start"].strftime(NAME_TIME_FORMAT)
    return (
        "Horizontal Grids:\n"
        "Name, H-Coord, nX, nY, dX, dY, X Min, Y Min\n"
        f"HGrid1, Lat-Long, {g.nx}, {g.ny}, {g.dx}, {g.dy}, {g.x_min}, {g.y_min}\n"
        "\n"
        "Vertical Grids:\n"
        "Name, Z-Coord, nZ, Z0, dZ\n"
        # One layer from the ground to `layer` m: centre layer/2, thickness `layer`.
        f"ZGrid1, m agl, 1, {layer / 2}, {layer}\n"
        "\n"
        "Temporal Grids:\n"
        "Name, nt, dt, t0\n"
        # A single output time at the chronologically earliest end of the window.
        f"TGridEndOfRun, 1, 01:00, {t0}\n"
    )


def domain_block(settings, timing):
    d = settings.domain
    return (
        "Domains:\n"
        "Name, H Unbounded?, H-Coord, X Min, X Max, Y Min, Y Max, Z Unbounded?, Z-Coord, Z Max, "
        "T Unbounded?, Start Time, End Time, Max Travel Time\n"
        f"Dispersion Domain, No, Lat-Long, {d.x_min}, {d.x_max}, {d.y_min}, {d.y_max}, No, m agl, "
        f"{d.z_max_m_agl}, Yes, , , {timing['max_travel_hours']}:00\n"
    )


def species_block():
    # Molecular weight 1 g/mol and a release of 1 g/s: mixing ratio output then converts to a
    # footprint without any molar-mass factor (footprints.footprint_from_mixing_ratio).
    return (
        "Species:\n"
        "Name, Category, Half Life, UV Loss Rate, Surface Resistance, Deposition Velocity, "
        "Molecular Weight, Material Unit\n"
        "INERT-TRACER, INERT, Stable, 0.00E+00, , 0.0, 1.0, g\n"
        "\n"
        "Species Uses:\n"
        "Species, On Particles?, On Fields?, Advect Field?\n"
        "INERT-TRACER, Yes, No, No\n"
    )


def dispersion_options_block(settings):
    return (
        "Sets of Dispersion Options:\n"
        "Max # Particles, Max # Full Particles, Velocity Memory Time, "
        "Mesoscale Velocity Memory Time, Inhomogeneous Time, Skew Time, Sync Time, DeltaOpt, "
        "Time of Fixed Met, Computational Domain, Deep Convection?, Radioactive Decay?, "
        "Agent Decay?, Dry Deposition?, Wet Deposition?, Turbulence?, Mesoscale Motions?, "
        "Chemistry?\n"
        f"{settings.max_particles}, 2, 00:00, 00:00, 00:00, 00:00, 00:{settings.sync_minutes:02d}, "
        "1, , Dispersion Domain, No, No, No, No, No, Yes, Yes, No\n"
    )


def particle_fate_block(settings):
    # One line per particle when it stops (left the domain, or reached the maximum travel time):
    # source number, travel time and position (dev/notes/name-behaviour.md); with fate_met, also
    # the met at that point (U/V/W Ambient ... Pressure ...).
    met = "Yes" if settings.fate_met else "No"
    return (
        "Output Requirements - Sets of Particle/Puff Information:\n"
        "Output Name, Particles?, Puffs?, Fate Info?, Met?, Mass?, Plume Rise?, "
        "Dispersion Scheme?, Puff Family?, H-Coord, Z-Coord, Sync?, Output Route\n"
        f"Particles, Yes, No, Yes, {met}, No, No, No, No, Lat-Long, {settings.fate_z_coord}, "
        "Yes, D\n"
    )


# NAME writes each output group's fields side by side, one column per field request, and stops
# with "too long an output line required for output group" if the line is too long - but only
# when it first writes, i.e. at the END of the run. 2,720 columns in one group failed after 8 h
# (stage 2 batch); 2,000 per group works (dev/notes/name-behaviour.md).
MAX_FIELD_COLUMNS_PER_GROUP = 2000
FIELDS_GROUP = "Fields_grid1"


def fields_output_group(index):
    """Output group of the release at 0-based position ``index`` in the plan.

    The first MAX_FIELD_COLUMNS_PER_GROUP releases go to "Fields_grid1" (the legacy name), the next
    to "Fields_grid1_2", and so on. NAME writes one fields file per group
    (``<group>_C1_T1_<time>.txt``); ``name_output.read_fields_folder`` reads them all.
    """
    group = index // MAX_FIELD_COLUMNS_PER_GROUP + 1
    return FIELDS_GROUP if group == 1 else f"{FIELDS_GROUP}_{group}"


def sources_and_requests_block(plan, settings, timing):
    """One cuboid source and one integrated mixing-ratio field request per release.

    Field requests are split into output groups of at most MAX_FIELD_COLUMNS_PER_GROUP releases
    (see ``fields_output_group``).
    """
    n_steps = timing["window_hours"] * 60 // settings.sync_minutes
    lines = [
        "Sources:",
        "Name, Start Time, Stop Time, # Particles, Source Strength, Shape, H-Coord, Z-Coord, "
        "Set of Locations, Location, dH-Metres?, dZ-Metres?, X, Y, Z, dX, dY, dZ, Angle, "
        "Top Hat, Plume Rise?, Temperature, Volume Flow Rate, Max Age, No Reflect?",
    ]
    for release in plan.itertuples():
        # Backwards run: the source "starts" at the later time and "stops" at the earlier one.
        lines.append(
            f"{release.source_name}, {release.release_start.strftime(NAME_TIME_FORMAT)}, "
            f"{release.release_end.strftime(NAME_TIME_FORMAT)}, "
            f"{release.particles_per_second}/sec, INERT-TRACER 1.0 g/s, Cuboid, Lat-Long, Pa, , , "
            f"No, No, {release.lon!r}, {release.lat!r}, {release.z_pa!r}, {release.dlon!r}, "
            f"{release.dlat!r}, {release.dz_pa!r}, 0.0, Yes, No, 0.0, 0.0, infinity, Yes"
        )
    lines += [
        "",
        "Output Requirements - Fields:",
        "Name, Quantity, Species, Source, H-Grid, Z-Grid, T-Grid, BL Average, T Av Or Int, "
        "Av Time, # Av Times, Sync?, Output Route, Across, Separate File, Output Format, "
        "Output Group",
    ]
    for index, release in enumerate(plan.itertuples()):
        number = release.source_name.removeprefix("SourceID_")
        lines.append(
            f"Req_H1Z1_{number}, Mixing Ratio, INERT-TRACER, {release.source_name}, HGrid1, "
            f"ZGrid1, TGridEndOfRun, No, Int, {timing['window_hours']}:00, {n_steps}, No, D, TZ, "
            f"T, IA, {fields_output_group(index)}"
        )
    return "\n".join(lines) + "\n"


def met_declarations_block(settings):
    """The NAME team's 8.7.1 met declarations for this met version, filled in for this run.

    The demo file has placeholders %MetRestoreScript%, %MetDir% and %TopogDir%. The minimum
    boundary-layer depth (first number in each "NWP Met Module Instances" row) is set to
    ``settings.min_bl_depth_m`` (legacy 40 m; demo 100 m).
    """
    path = DEMO_UTILITIES / "MetDeclarations" / settings.met_version.declaration_file
    text = path.read_text()
    text = (
        text.replace("%MetRestoreScript%", settings.met_restore_script)
        .replace("%MetDir%", settings.met_folder)
        .replace("%TopogDir%", settings.topography_folder)
    )
    if "%" in text:
        raise ValueError(f"Unfilled placeholder left in {path.name}")

    lines = text.splitlines()
    in_met_instances = False
    for i, line in enumerate(lines):
        if line.startswith("NWP Met Module Instances:"):
            in_met_instances = True
            continue
        if in_met_instances:
            if not line.strip():
                in_met_instances = False
            elif not line.startswith("Name,"):
                cells = line.split(",")
                cells[1] = f" {settings.min_bl_depth_m}"
                lines[i] = ",".join(cells)
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------------------------
# The whole file
# ---------------------------------------------------------------------------------------------


def name_input_text(plan, settings):
    """Text of the NAME input file for one day's release plan."""
    timing = run_timing(plan, settings)
    blocks = [
        "! NAME input file written by name_satellite (https://github.com/mrghg/name-satellite)\n",
        main_options_block(settings, len(plan)),
        restart_and_cases_block(),
        openmp_block(settings),
        timer_block(settings),
        files_and_coordinates_block(settings),
        grids_block(settings, timing),
        domain_block(settings, timing),
        species_block(),
        dispersion_options_block(settings),
        particle_fate_block(settings),
        sources_and_requests_block(plan, settings, timing),
        met_declarations_block(settings),
    ]
    if settings.release_met:
        soundings = release_met_soundings(plan)
        blocks.append(release_met_blocks(soundings, settings.grid.dx, settings.grid.dy))
    return "\n".join(blocks)


def write_name_input(path, plan, settings):
    """Write the NAME input file; returns the path."""
    path = Path(path)
    path.write_text(name_input_text(plan, settings))
    return path
