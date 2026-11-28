"""Prepare a month of future-method NAME runs (and profiling variants) from Level-2 data.

Where this sits (docs/future/scoping.md)
----------------------------------------
    L2 files  ->  **future.setup**: per day, observation file + surface-pressure run + main job
    ->  on LOTUS2: met restore -> surface-pressure run -> main job (writes the NAME input with
    ``main-input``, runs NAME with ``future.sampler`` beside it)  ->  ``future.processing``

Steps 1-7 are those of the main workflow (``prepare``): soundings, the model surface pressure
(D19), modelled levels (D21), particle rates (D22), release groups (D20). What differs:

* **No side boundaries**: particles travel for ``MAX_AGE_HOURS`` (240 h, the proposed default)
  wherever they go; the top stays at 20 km agl (D21). Their end positions (the patched NAME's fate
  records, X20) give the boundary term.
* **Nested output grids** (``future.grids``): per release group a 0.1° box (±5°) and a 0.5° box
  (±25°), snapped to the 2.5° lines; one global 2.5° grid for every release. One output group per
  grid (NAME needs consistent grids within a group).
* **Release groups** (D20) use 0.1° cells of the global grid instead of the South America grid.
* **Physics as NAME recommends** (NAME documentation A01 §1.7.5-1.7.6):
  - deep convection ``New``;
  - the near-source turbulence scheme (velocity memory, height-dependent turbulence profiles) for
    each particle's first ``NEAR_SOURCE_HOURS`` of travel, then the long-range scheme. Skewed
    turbulence off (``SKEW_HOURS`` = 0). ``MAX_FULL_PARTICLES`` must cover every particle the run
    releases, not just those still within that time of release: NAME stopped ("array is full") on
    days releasing 1.36-1.44 M particles with a limit of 1 M, although at most 0.37 M were within an
    hour of release (Jan 2013), so it is set to ``MAX_PARTICLES``;
  - boundary-layer depth limits as the NAME team's met declarations (minimum 100 m; the main
    workflow keeps legacy's 40 m).
* No release-point met requests (not used by ``future.processing``).
* **Layouts** (``Options.layout``): how NAME's output maps to soundings.
  - ``species`` (default): one species per sounding. Each release emits the species of every
    sounding sharing its release group, at that sounding's weight for the release's level
    (pressure weight x averaging kernel, g/s); three field requests per sounding (by species).
    NAME then outputs each sounding's column directly. At most ``MAX_SPECIES_PER_RUN`` soundings
    per run.
  - ``sources``: one species; three field requests per release, each restricted to it; the
    levels are weighted in ``future.processing``. Gives each level's footprint, but run time
    grows with particles x field requests in the run, so it is much slower for large runs.
  Both give the same column footprints (to NAME's output precision; 7 Jan 2013, 42 soundings:
  cells within 3e-6 of each field's maximum). See dev/notes/name-behaviour.md and
  docs/future/run-sizing.md for the timing.
* Soundings whose model surface pressure is outside 0.9-1.1 × the retrieval's (steep terrain: the
  model ground more than ~1 km off) are dropped and listed in ``dropped_soundings.csv`` (D27)
  instead of stopping the day, as the main workflow does.

Profiling: every NAME job runs ``future.sampler`` beside NAME (CPU per thread, bytes read and
written, memory, met progress every 5 s, in ``profile.csv``), and NAME's timer summary is on.
With ``--profiling``, extra runs of days already in the month change one thing each
(``PROFILING_VARIANTS``): threads, prefetch, detailed timers, convection, the near-source scheme,
and three days in one run. They reuse their days' observation files and surface-pressure runs.
``split`` makes the opposite test: one day's release groups shared between N runs (``split_plan``),
as run time grows with particles x releases in a run.

Command line
------------
::

    uv run python -m name_satellite.future.setup PRODUCT SITE YYYY-MM WORK [--test] [--profiling]
        [--days=D1,D2,...]
    uv run python -m name_satellite.future.setup main-input RUN_FOLDER
    uv run python -m name_satellite.future.setup split MONTH_FOLDER YYYY-MM-DD[,...] N
        [--parts=1,2] [--threads=T] [--output-grids=coarse] [--label=NAME]

The month folder is ``WORK/<SITE>_<PRODUCT>_<YYYYMM>_future/``; submit with its ``submit.sh``.
``--test``: 1 hour of travel, debug queue. ``main-input`` is what each main job runs first.
``split`` writes ``profile_<label>_part<k>`` folders (default label ``split<N>_<YYYYMMDD>``) for
already-prepared days, each with its job (submit with ``sbatch``); options as listed above.

Files: as ``prepare`` (``summary.csv``, ``restore_met.sbatch``, ``submit.sh``; per day
``obs_*.nc``, ``day.json``, ``surface_pressure/``, ``job.sbatch``; written on LOTUS2:
``surface_pressure.csv``, ``release_groups.csv``, ``release_plan.csv``, ``input.txt``,
``output/``), plus ``grids.csv`` (each release group's boxes) and ``profile.csv``. A profiling
variant folder (``profile_<name>``) has ``day.json`` naming the day folders it uses.
"""

import dataclasses
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from name_satellite import name_input
from name_satellite.batching import combine_release_plans
from name_satellite.domains import SITES
from name_satellite.future.grids import (
    COARSE,
    FINE,
    MEDIUM,
    Box,
    box_shape,
    global_box,
    name_grid_line,
    snapped_box,
)
from name_satellite.l2_files import L2FileMissing, download_l2
from name_satellite.lotus import REPO, write_job_script, write_met_restore_job
from name_satellite.met import PATCHED_NAME_DIR, met_version_for_period
from name_satellite.name_input import NAME_TIME_FORMAT, OutputGrid, RunSettings, run_timing
from name_satellite.observations import read_l2, to_obs, write_obs_file
from name_satellite.prepare import (
    L2_CACHE,
    MET_ROOT,
    grouped_release_plan,
    release_boxes,
    select_soundings,
)
from name_satellite.release_met import MAX_SOUNDINGS, release_met_soundings
from name_satellite.releases import (
    MODELLED_LEVELS,
    PARTICLE_RATES,
    build_release_plan,
    check_modelled_levels,
    scale_to_model_surface,
)
from name_satellite.surface_pressure import (
    implausible_surface_pressure,
    read_surface_pressure,
    surface_pressure_input_text,
    surface_pressure_settings,
    surface_pressure_table,
)

MAX_AGE_HOURS = 240
SYNC_MINUTES = 6
# 2 threads: met reading does not speed up with threads, and a ~10.5 GB run on 2 cores matches
# LOTUS2's 8 GB per core, so 2 threads use the least of a node per run (docs/future/run-sizing.md).
THREADS = 2
# Memory to request for the main NAME job, by met version. Mk6/Mk7 runs peaked at 10.3-10.8 GB.
# Newer versions have finer grids (Mk8/9: 2.2x the points, Mk10-12: 6x) and are not yet measured:
# 48 GB until they are (docs/future/run-sizing.md, "Memory and met version").
JOB_MEMORY = {"UMG_Mk6PT": "16G", "UMG_Mk7PT": "16G"}
UNMEASURED_JOB_MEMORY = "48G"
NEAR_SOURCE_HOURS = 1.0  # velocity memory + inhomogeneous turbulence (A01 §1.7.5.5: 30 min-1 h)
SKEW_HOURS = 0.0
DEEP_CONVECTION = "New"
MIN_BL_DEPTH_M = 100.0  # the NAME team's met declarations
MAX_PARTICLES = 7_000_000
MAX_FULL_PARTICLES = MAX_PARTICLES  # see the module docstring
MAX_FIELD_REQUIREMENTS = 20_000
MAX_FIELD_OUTPUT_GROUPS = 1000
# NAME stops at once ("... is not large enough") with more than 3,000 field requests of this kind
# in a run, whatever their horizontal or time grids or output groups (tests 26 Nov 2026,
# dev/scripts/name_tests/field_request_limit.py). Three per release: at most 1,000 releases.
MAX_FIELD_REQUESTS_PER_RUN = 3000
REQUESTS_PER_RELEASE = {"nested": 3, "coarse": 1}  # per Options.output_grids
# Species particles can carry in NAME 8.7.1 (not yet reached in a run): in the "species" layout,
# at most this many soundings per run.
MAX_SPECIES_PER_RUN = 216
SAMPLER = REPO / "src" / "name_satellite" / "future" / "sampler.py"

# Release groups share releases within a 0.1° cell of the global grid (D20, with this grid).
GLOBAL_FINE_GRID = OutputGrid(nx=3600, ny=1800, dx=0.1, dy=0.1, x_min=-179.95, y_min=-89.95)


@dataclasses.dataclass(frozen=True)
class Options:
    """What a run may vary (defaults: the future method's settings)."""

    threads: int = THREADS
    deep_convection: str = DEEP_CONVECTION
    near_source_hours: float = NEAR_SOURCE_HOURS
    skew_hours: float = SKEW_HOURS
    prefetch: bool = True
    detailed_timers: bool = False
    test: bool = False
    # "nested": fine, medium and coarse requests per release; "coarse": the global 2.5° grid only
    # (a profiling test of how run time depends on field requests; not processable)
    output_grids: str = "nested"
    # How releases map to footprints in NAME (see the module docstring, "Layouts"):
    # "species": one species per sounding, its weights as the releases' strengths (default);
    # "sources": one species, field requests per release, weights applied in processing.
    layout: str = "species"


# Profiling variants: name -> (days used, options). Each changes one thing from the defaults.
PROFILING_VARIANTS = {
    "threads8": (["2013-01-10"], {"threads": 8}),
    "noprefetch": (["2013-01-11"], {"prefetch": False}),
    "detailedtimers": (["2013-01-12"], {"detailed_timers": True}),
    "noconvection": (["2013-01-13"], {"deep_convection": "No"}),
    "nonearsource": (["2013-01-14"], {"near_source_hours": 0.0}),
    "batch3days": (["2013-01-03", "2013-01-04", "2013-01-05"], {}),  # 528 releases
}


def future_site(site):
    """The site with release groups on 0.1° cells of the global grid (D20 with that grid)."""
    domain = dataclasses.replace(site.domain, grid=GLOBAL_FINE_GRID)
    return dataclasses.replace(site, domain=domain, release_cell="grid")


def hhmm(hours):
    """Hours as NAME's HH:MM."""
    minutes = round(hours * 60)
    return f"{minutes // 60:02d}:{minutes % 60:02d}"


def run_settings(site, run_folder, first_day, last_day, options):
    """NAME settings for a run whose releases fall on ``first_day``..``last_day``.

    The met version must cover from MAX_AGE_HOURS (+1 day) before the first day to the day after
    the last (``met.met_version_for_period``, D25).
    """
    hours = 1 if options.test else MAX_AGE_HOURS
    met_version = met_version_for_period(
        first_day - pd.Timedelta(hours=hours + 24), last_day + pd.Timedelta(days=1)
    )
    return RunSettings(
        run_name=run_folder.name,
        output_folder=str(run_folder / "output"),
        grid=GLOBAL_FINE_GRID,
        domain=site.domain.dispersion,  # only its top is used (unbounded horizontally)
        met_version=met_version,
        met_folder=str(MET_ROOT / met_version.name),
        max_age_hours=hours,
        sync_minutes=SYNC_MINUTES,
        max_particles=MAX_PARTICLES,
        max_field_requirements=MAX_FIELD_REQUIREMENTS,
        max_field_output_groups=MAX_FIELD_OUTPUT_GROUPS,
        threads=options.threads,
        min_bl_depth_m=MIN_BL_DEPTH_M,
        release_met=False,
        use_timers=True,
    )


# ---------------------------------------------------------------------------------------------
# NAME input blocks that differ from the main workflow's (name_input)
# ---------------------------------------------------------------------------------------------


def group_boxes(plan):
    """One row per release group: its fine and medium boxes (degrees) and NAME grid names.

    Groups are numbered 1, 2, ... in plan order (first appearance).
    """
    rows = []
    for number, group in enumerate(dict.fromkeys(plan["release_group"]), 1):
        first = plan[plan["release_group"] == group].iloc[0]
        fine = snapped_box(first["lon"], first["lat"], FINE.half_width_deg)
        medium = snapped_box(first["lon"], first["lat"], MEDIUM.half_width_deg)
        rows.append(
            {"release_group": group, "group_number": number, "fine_grid": f"HGridF{number}",
             "medium_grid": f"HGridM{number}", **{f"fine_{k}": v for k, v in vars(fine).items()},
             **{f"medium_{k}": v for k, v in vars(medium).items()}}
        )  # fmt: skip
    return pd.DataFrame(rows)


def box_of(row, level):
    """The ``level`` ("fine"/"medium") box of a ``group_boxes`` row (dict or Series)."""
    return Box(row[f"{level}_west"], row[f"{level}_east"], row[f"{level}_south"],
               row[f"{level}_north"])  # fmt: skip


def timer_block(options):
    summary = "No" if options.detailed_timers else "Yes"
    return f"Timer Options:\nUse Timers?, Summary Only?\nYes, {summary}\n"


def grids_block(boxes, settings, timing, output_grids="nested"):
    """Horizontal grids (every group's fine and medium box unless ``output_grids`` is "coarse",
    the global coarse grid), and the vertical and temporal grids of the main workflow."""
    lines = ["Horizontal Grids:", "Name, H-Coord, nX, nY, dX, dY, X Min, Y Min"]
    for row in boxes.itertuples() if output_grids == "nested" else []:
        lines.append(name_grid_line(row.fine_grid, box_of(row._asdict(), "fine"), FINE.cell_deg))
        lines.append(
            name_grid_line(row.medium_grid, box_of(row._asdict(), "medium"), MEDIUM.cell_deg)
        )
    lines.append(name_grid_line("HGridC", global_box(), COARSE.cell_deg))
    standard = name_input.grids_block(settings, timing)
    return "\n".join(lines) + "\n\n" + standard[standard.index("Vertical Grids:") :]


def domain_block(settings, timing):
    """No side boundaries (NAME refuses an H-Coord or bounds then); top as the main workflow."""
    return (
        "Domains:\n"
        "Name, H Unbounded?, H-Coord, X Min, X Max, Y Min, Y Max, Z Unbounded?, Z-Coord, Z Max, "
        "T Unbounded?, Start Time, End Time, Max Travel Time\n"
        f"Dispersion Domain, Yes, , , , , , No, m agl, {settings.domain.z_max_m_agl}, Yes, , , "
        f"{timing['max_travel_hours']}:00\n"
    )


def dispersion_options_block(settings, options):
    """Particle limits, turbulence scheme, sync time, convection (NAME docs B01, A01 §1.7)."""
    if options.deep_convection not in ("No", "Old", "New"):
        raise ValueError(f"Deep convection must be No, Old or New, not {options.deep_convection}")
    if not 0 <= options.skew_hours <= options.near_source_hours:
        raise ValueError("Need 0 <= skew time <= velocity memory time (NAME docs B01)")
    near = hhmm(options.near_source_hours)
    full = MAX_FULL_PARTICLES if options.near_source_hours > 0 else 2  # NAME docs A01 §1.5.6
    return (
        "Sets of Dispersion Options:\n"
        "Max # Particles, Max # Full Particles, Velocity Memory Time, "
        "Mesoscale Velocity Memory Time, Inhomogeneous Time, Skew Time, Sync Time, DeltaOpt, "
        "Time of Fixed Met, Computational Domain, Deep Convection?, Radioactive Decay?, "
        "Agent Decay?, Dry Deposition?, Wet Deposition?, Turbulence?, Mesoscale Motions?, "
        "Chemistry?\n"
        f"{settings.max_particles}, {full}, {near}, 00:00, {near}, {hhmm(options.skew_hours)}, "
        f"00:{settings.sync_minutes:02d}, 1, , Dispersion Domain, {options.deep_convection}, "
        "No, No, No, No, Yes, Yes, No\n"
    )


def sources_and_requests_block(plan, boxes, settings, timing, output_grids="nested"):
    """Sources as the main workflow; per source three field requests (its group's fine and medium
    grids, the global coarse grid; only the last if ``output_grids`` is "coarse"), one output group
    per grid."""
    standard = name_input.sources_and_requests_block(plan, settings, timing)
    sources = standard[: standard.index("Output Requirements - Fields:")]
    n_steps = timing["window_hours"] * 60 // settings.sync_minutes
    by_group = boxes.set_index("release_group")
    lines = [
        "Output Requirements - Fields:",
        "Name, Quantity, Species, Source, H-Grid, Z-Grid, T-Grid, BL Average, T Av Or Int, "
        "Av Time, # Av Times, Sync?, Output Route, Across, Separate File, Output Format, "
        "Output Group",
    ]
    for release in plan.itertuples():
        number = release.source_name.removeprefix("SourceID_")
        group = by_group.loc[release.release_group]
        grids = (group["fine_grid"], group["medium_grid"], "HGridC")
        for grid in grids if output_grids == "nested" else grids[-1:]:
            lines.append(
                f"Req_{grid}_{number}, Mixing Ratio, INERT-TRACER, {release.source_name}, {grid}, "
                f"ZGrid1, TGridEndOfRun, No, Int, {timing['window_hours']}:00, {n_steps}, No, D, "
                f"TZ, T, IA, Fields_{grid}"
            )
    return sources + "\n".join(lines) + "\n"


def species_names(weights):
    """NAME species name of each sounding (``SND<n>``, in order of first appearance): a Series
    indexed by sounding_id."""
    soundings = list(dict.fromkeys(weights["sounding_id"]))
    return pd.Series([f"SND{k}" for k in range(1, len(soundings) + 1)],
                     index=pd.Index(soundings, name="sounding_id"), name="species")  # fmt: skip


def species_main_options_block(settings, n_sources, n_species):
    """``Main Options`` as ``name_input.main_options_block``, with room for every species."""
    return (
        "Main Options:\n"
        "Absolute or Relative Time?, Fixed Met?, Flat Earth?, Run Name, Random Seed, "
        "Max # Species, Max # Sources, Max # Field Reqs, Max # Field Output Groups, Backwards?\n"
        f"Absolute, No, No, {settings.run_name}, {settings.random_seed}, {n_species}, "
        f"{n_sources}, {settings.max_field_requirements}, {settings.max_field_output_groups}, "
        "Yes\n"
    )


def species_block(species):
    """One inert tracer per sounding, as ``name_input.species_block``'s INERT-TRACER (molecular
    weight 1 g/mol, so a release of w g/s is w mol/s)."""
    lines = ["Species:",
             "Name, Category, Half Life, UV Loss Rate, Surface Resistance, Deposition Velocity, "
             "Molecular Weight, Material Unit"]  # fmt: skip
    lines += [f"{name}, INERT, Stable, 0.00E+00, , 0.0, 1.0, g" for name in species]
    lines += ["", "Species Uses:", "Species, On Particles?, On Fields?, Advect Field?"]
    lines += [f"{name}, Yes, No, No" for name in species]
    return "\n".join(lines) + "\n"


def species_sources_and_requests_block(plan, boxes, weights, species, settings, timing,
                                       output_grids="nested"):  # fmt: skip
    """Sources and field requests of the "species" layout.

    Each release emits the species of every sounding sharing its release group, at that
    sounding's weight for the release's level (g/s, i.e. mol/s: molecular weight 1). A sounding's
    field is then already Σ_j weight_j x (level j's field at 1 g/s), the column footprint once
    converted (``footprints.footprint_from_mixing_ratio``). Three field requests per sounding (its
    group's fine and medium grids and the global coarse grid; only the last if ``output_grids``
    is "coarse"), by species, not restricted to a source; one output group per grid.

    The source lines are those of ``name_input.sources_and_requests_block`` except for the
    strengths. Weights are written to 7 significant digits.
    """
    per_group = weights.groupby("release_group", sort=False)["sounding_id"].unique()
    n_columns = max(len(v) for v in per_group)
    weight = weights.set_index(["sounding_id", "level_number"])["weight"]
    lines = [
        "Sources:",
        "Name, Start Time, Stop Time, # Particles, "
        + ", ".join(["Source Strength"] * n_columns)
        + ", Shape, H-Coord, Z-Coord, Set of Locations, Location, dH-Metres?, dZ-Metres?, X, Y, "
        "Z, dX, dY, dZ, Angle, Top Hat, Plume Rise?, Temperature, Volume Flow Rate, Max Age, "
        "No Reflect?",
    ]
    for release in plan.itertuples():
        strengths = [f"{species[s]} {weight[(s, release.level_number)]:.6e} g/s"
                     for s in per_group[release.release_group]]  # fmt: skip
        strengths += [""] * (n_columns - len(strengths))
        # Backwards run: the source "starts" at the later time and "stops" at the earlier one.
        lines.append(
            f"{release.source_name}, {release.release_start.strftime(NAME_TIME_FORMAT)}, "
            f"{release.release_end.strftime(NAME_TIME_FORMAT)}, "
            f"{release.particles_per_second}/sec, {', '.join(strengths)}, Cuboid, Lat-Long, Pa, "
            f", , No, No, {release.lon!r}, {release.lat!r}, {release.z_pa!r}, {release.dlon!r}, "
            f"{release.dlat!r}, {release.dz_pa!r}, 0.0, Yes, No, 0.0, 0.0, infinity, Yes"
        )
    n_steps = timing["window_hours"] * 60 // settings.sync_minutes
    lines += [
        "",
        "Output Requirements - Fields:",
        "Name, Quantity, Species, Source, H-Grid, Z-Grid, T-Grid, BL Average, T Av Or Int, "
        "Av Time, # Av Times, Sync?, Output Route, Across, Separate File, Output Format, "
        "Output Group",
    ]
    by_group = boxes.set_index("release_group")
    group_of = weights.drop_duplicates("sounding_id").set_index("sounding_id")["release_group"]
    for sounding, name in species.items():
        group = by_group.loc[group_of[sounding]]
        grids = (group["fine_grid"], group["medium_grid"], "HGridC")
        for grid in grids if output_grids == "nested" else grids[-1:]:
            lines.append(
                f"Req_{grid}_{name}, Mixing Ratio, {name}, , {grid}, ZGrid1, TGridEndOfRun, No, "
                f"Int, {timing['window_hours']}:00, {n_steps}, No, D, TZ, T, IA, Fields_{grid}"
            )
    return "\n".join(lines) + "\n"


def met_declarations_block(settings, options):
    """The NAME team's declarations (as ``name_input``), with ``Prefetch? No`` if asked."""
    text = name_input.met_declarations_block(settings)
    if options.prefetch:
        return text
    lines = text.splitlines()
    in_instances = False
    for i, line in enumerate(lines):
        if line.startswith("NWP Met Module Instances:"):
            in_instances = True
            continue
        if in_instances:
            if not line.strip():
                in_instances = False
            elif line.startswith("Name,"):
                lines[i] = line.rstrip().rstrip(",") + ", Prefetch?"
            else:
                lines[i] = line.rstrip().rstrip(",") + ", No"
    return "\n".join(lines) + "\n"


def input_text(plan, boxes, settings, options, weights=None):
    """The whole NAME input file of a future-method run.

    weights: the soundings' weights (``weights_table``) of the plan's release groups; needed for
    the "species" layout (one species per sounding), unused for "sources".
    """
    if options.layout == "species":
        if weights is None:
            raise ValueError("The species layout needs the soundings' weights")
        n_species = weights["sounding_id"].nunique()
        if n_species > MAX_SPECIES_PER_RUN:
            raise ValueError(f"{n_species} soundings in one run; particles can carry at most "
                             f"{MAX_SPECIES_PER_RUN} species")  # fmt: skip
        n_requests = REQUESTS_PER_RELEASE[options.output_grids] * n_species
    elif options.layout == "sources":
        n_requests = REQUESTS_PER_RELEASE[options.output_grids] * len(plan)
    else:
        raise ValueError(f"Unknown layout {options.layout!r}")
    if n_requests > MAX_FIELD_REQUESTS_PER_RUN:
        raise ValueError(f"The run needs {n_requests} field "
                         f"requests; NAME allows {MAX_FIELD_REQUESTS_PER_RUN} per run")  # fmt: skip
    for column in ("fine", "medium"):
        too_big = [box_shape(box_of(r, column), (FINE if column == "fine" else MEDIUM).cell_deg)
                   for _, r in boxes.iterrows()]  # fmt: skip
        limit = (FINE if column == "fine" else MEDIUM).block_cells
        if max(max(s) for s in too_big) > limit:
            raise ValueError(f"A {column} box is larger than its {limit}-cell block")
    timing = run_timing(plan, settings)
    if options.layout == "species":
        species = species_names(weights)
        main_options = species_main_options_block(settings, len(plan), len(species))
        species_text = species_block(species)
        sources_and_requests = species_sources_and_requests_block(
            plan, boxes, weights, species, settings, timing, options.output_grids
        )
    else:
        main_options = name_input.main_options_block(settings, len(plan))
        species_text = name_input.species_block()
        sources_and_requests = sources_and_requests_block(plan, boxes, settings, timing,
                                                          options.output_grids)  # fmt: skip
    blocks = [
        "! NAME input file: future method (name_satellite.future.setup)\n",
        main_options,
        name_input.restart_and_cases_block(),
        name_input.openmp_block(settings),
        timer_block(options),
        name_input.files_and_coordinates_block(settings),
        grids_block(boxes, settings, timing, options.output_grids),
        domain_block(settings, timing),
        species_text,
        dispersion_options_block(settings, options),
        name_input.particle_fate_block(settings),
        sources_and_requests,
        met_declarations_block(settings, options),
    ]
    return "\n".join(blocks)


# ---------------------------------------------------------------------------------------------
# Days, variants and the month
# ---------------------------------------------------------------------------------------------


def column_weights(obs, n_levels):
    """Pressure weight x averaging kernel of each sounding's modelled levels: array (sounding,
    level), dimensionless. Each sounding's column footprint is Σ_j weight_j x level j's footprint.

    Raises:
        ValueError if a weight is missing or negative.
    """
    weights = (obs["pressure_weights"].values[:, :n_levels]
               * obs["xch4_averaging_kernel"].values[:, :n_levels]).astype(float)  # fmt: skip
    if not np.isfinite(weights).all():
        raise ValueError("Missing pressure weights or averaging kernels in the modelled levels")
    if (weights < 0).any():
        raise ValueError("Negative pressure weight x averaging kernel in the modelled levels")
    return weights


def weights_table(obs, members, n_levels):
    """Long table of the weights: sounding_id, release_group, level_number, weight."""
    weights = column_weights(obs, n_levels)
    group_of = members.set_index("sounding_id")["release_group"]
    rows = []
    for s, sounding in enumerate(str(i) for i in obs["sounding_id"].values):
        for level in range(1, n_levels + 1):
            rows.append((sounding, group_of[sounding], level, weights[s, level - 1]))
    return pd.DataFrame(rows, columns=["sounding_id", "release_group", "level_number", "weight"])


def final_plan(day_folder, site, product, write_tables=True):
    """A day's release plan with levels scaled to the model surface (D19): as
    ``prepare.write_main_input``, with release groups on the global 0.1° grid. Returns the plan
    and the soundings' weights (``weights_table``).

    With ``write_tables``, writes ``surface_pressure.csv`` and ``release_groups.csv`` in the day
    folder (only the day's own job does; a profiling variant must not write into a day folder
    that the day's job may be writing at the same time).
    """
    day_folder = Path(day_folder)
    day = pd.Timestamp(json.loads((day_folder / "day.json").read_text())["day"])
    with xr.open_dataset(day_folder / f"obs_{site.name}_{day:%Y%m%d}.nc") as f:
        obs = f.load()
    folder = day_folder / "surface_pressure"
    soundings = pd.read_csv(folder / "soundings.csv", parse_dates=["time"],
                            dtype={"sounding_id": str})  # fmt: skip
    surface = surface_pressure_table(obs, read_surface_pressure(folder / "output", soundings),
                                     check_ratio=False)  # fmt: skip
    # Soundings whose model ground is more than ~10% off the retrieval's surface pressure (steep
    # terrain) are dropped and listed (D27); future.processing skips exactly those.
    bad = implausible_surface_pressure(surface).to_numpy()
    if write_tables:
        surface.to_csv(day_folder / "surface_pressure.csv", index=False)
        surface[bad].to_csv(day_folder / "dropped_soundings.csv", index=False)
    if bad.any():
        print(f"{day_folder}: dropped {bad.sum()} soundings (implausible surface pressure, D27)")
    obs = obs.isel(sounding_id=np.flatnonzero(~bad))
    surface = surface[~bad].reset_index(drop=True)
    levels = obs["pressure_levels"].values
    scaled = np.array([scale_to_model_surface(levels[p], surface["model_hpa"].iloc[p])
                       for p in range(len(surface))])  # fmt: skip
    n_levels = MODELLED_LEVELS[product]
    check_modelled_levels(obs, scaled, surface["domain_top_hpa"].to_numpy(), n_levels)
    plan, members = grouped_release_plan(obs, future_site(site), n_levels,
                                         pressure_levels_hpa=scaled,
                                         particle_rates=PARTICLE_RATES[product])  # fmt: skip
    if write_tables:
        members.to_csv(day_folder / "release_groups.csv", index=False)
    return plan, weights_table(obs, members, n_levels)


def write_main_input(run_folder):
    """Write a run's release plan, grids and NAME input (what each main job runs first).

    A day folder plans its own day; a profiling variant combines the plans of the day folders its
    ``day.json`` names (``batching.combine_release_plans``).
    """
    run_folder = Path(run_folder)
    info = json.loads((run_folder / "day.json").read_text())
    site = SITES[info["site"]]
    options = Options(**info["options"])
    day_folders = [Path(p) for p in info["day_folders"]]
    planned = {f.name: final_plan(f, site, info["product"], write_tables=(f == run_folder))
               for f in day_folders}  # fmt: skip
    plans = {name: plan for name, (plan, _) in planned.items()}
    plan = plans[day_folders[0].name] if len(plans) == 1 else combine_release_plans(plans)
    if "part" in info:
        plan = split_plan(plan, *info["part"])
    weights = pd.concat([w for _, w in planned.values()], ignore_index=True)
    weights = weights[weights["release_group"].isin(set(plan["release_group"]))]
    boxes = group_boxes(plan)
    plan.to_csv(run_folder / "release_plan.csv", index=False)
    boxes.to_csv(run_folder / "grids.csv", index=False)
    if options.layout == "species":
        species = species_names(weights)
        species.to_csv(run_folder / "species.csv")
    days = [pd.Timestamp(json.loads((f / "day.json").read_text())["day"]) for f in day_folders]
    settings = run_settings(site, run_folder, min(days), max(days), options)
    (run_folder / "input.txt").write_text(input_text(plan, boxes, settings, options, weights))
    print(f"{run_folder}: {len(boxes)} release groups, {len(plan)} releases")
    return plan


def add_sampler(job_path, run_folder, input_path):
    """Run ``future.sampler`` beside NAME in the job (with the jaspy python3: stdlib only)."""
    job = job_path.read_text()
    run_line = next(line for line in job.splitlines() if line.startswith("time mpirun"))
    start = (f"python3 -I {SAMPLER} {input_path} {run_folder / 'profile.csv'} &\n"
             "sampler_pid=$!\n")  # fmt: skip
    if job.count(run_line) != 1 or job.count("status=$?") != 1:
        raise ValueError(f"Unexpected job script layout in {job_path}")
    job = job.replace(run_line, start + run_line)
    job = job.replace("status=$?", "status=$?\nwait $sampler_pid", 1)
    job_path.write_text(job)


def write_main_job(run_folder, job_name, options, met_version):
    """The main job: ``main-input``, then the patched NAME with the sampler beside it.

    met_version: the run's ``met.MetVersion`` (sets the memory request, ``JOB_MEMORY``).
    """
    memory = JOB_MEMORY.get(met_version.name, UNMEASURED_JOB_MEMORY)
    write = (f"cd {REPO}\n~/.local/bin/uv run --no-sync python -m name_satellite.future.setup "
             f"main-input {run_folder} || exit 1")  # fmt: skip
    qos, time = ("debug", "00:30:00") if options.test else ("high", "24:00:00")
    job = write_job_script(
        run_folder / "job.sbatch", run_folder, run_folder / "input.txt",
        threads=options.threads, qos=qos, time=time, memory=memory, job_name=job_name,
        before_name=write, name_dir=PATCHED_NAME_DIR,
    )  # fmt: skip
    add_sampler(job, run_folder, run_folder / "input.txt")


def prepare_day(product, site, day, month_folder, options):
    """Observation file, surface-pressure run and main job of one day (as ``prepare``)."""
    row = {"day": f"{day:%Y-%m-%d}", "status": "", "l2_soundings": 0, "soundings": 0, "sources": 0}
    try:
        l2_path = download_l2(product, day, L2_CACHE)
    except L2FileMissing:
        row["status"] = "no L2 file"
        return row
    table = read_l2(l2_path)
    row["l2_soundings"] = table.sizes["sounding"]
    table = select_soundings(table, site)
    row["soundings"] = table.sizes["sounding"]
    if row["soundings"] == 0:
        row["status"] = "no soundings in the box"
        return row
    if row["soundings"] > MAX_SOUNDINGS:
        raise ValueError(f"{day:%Y-%m-%d}: {row['soundings']} soundings > {MAX_SOUNDINGS}")

    run_folder = month_folder / f"{day:%Y%m%d}"
    (run_folder / "output").mkdir(parents=True)
    write_obs_file(table, run_folder / f"obs_{site.name}_{day:%Y%m%d}.nc")
    info = {"product": product, "site": site.name, "day": f"{day:%Y-%m-%d}",
            "options": dataclasses.asdict(options), "day_folders": [str(run_folder)]}  # fmt: skip
    (run_folder / "day.json").write_text(json.dumps(info, indent=2) + "\n")

    obs = to_obs(table)
    n_levels = MODELLED_LEVELS[product]
    sounding_plan = build_release_plan(obs, n_levels=n_levels,
                                       footprint_diameter_km=site.footprint_km,
                                       box_degrees=release_boxes(obs))  # fmt: skip
    plan, members = grouped_release_plan(obs, future_site(site), n_levels,
                                         particle_rates=PARTICLE_RATES[product])  # fmt: skip
    settings = run_settings(site, run_folder, day, day, options)

    # The met-only surface-pressure run (stock NAME), as the main workflow.
    folder = run_folder / "surface_pressure"
    (folder / "output").mkdir(parents=True)
    soundings = release_met_soundings(sounding_plan)
    soundings.to_csv(folder / "soundings.csv", index=False)
    sp_settings = surface_pressure_settings(dataclasses.replace(settings, max_age_hours=1),
                                            folder / "output")  # fmt: skip
    (folder / "input.txt").write_text(surface_pressure_input_text(soundings, sp_settings))
    write_job_script(folder / "job.sbatch", folder, folder / "input.txt", threads=1,
                     qos="debug" if options.test else "standard", time="00:30:00", memory="8G",
                     job_name=f"fut_{day:%Y%m%d}_ps")  # fmt: skip
    write_main_job(run_folder, f"fut_{day:%Y%m%d}", options, settings.met_version)

    timing = run_timing(plan, settings)
    row.update(status="prepared", release_groups=members["release_group"].nunique(),
               sources=len(plan), met_version=settings.met_version.name,
               met_start=f"{timing['window_start']:%Y-%m-%dT%H:%M}",
               met_end=f"{plan['release_start'].max():%Y-%m-%dT%H:%M}")  # fmt: skip
    return row


def prepare_variant(name, days, changes, month_folder, product, site, base_options):
    """A profiling run: the named days (already prepared) with one option changed."""
    folder = month_folder / f"profile_{name}"
    (folder / "output").mkdir(parents=True)
    options = dataclasses.replace(base_options, **changes)
    day_folders = [month_folder / f"{pd.Timestamp(d):%Y%m%d}" for d in days]
    info = {"product": product, "site": site.name, "day": days[0],
            "options": dataclasses.asdict(options), "day_folders": [str(f) for f in day_folders],
            "variant": name, "changes": changes}  # fmt: skip
    (folder / "day.json").write_text(json.dumps(info, indent=2) + "\n")
    dates = pd.to_datetime(days)
    settings = run_settings(site, folder, dates.min(), dates.max(), options)
    write_main_job(folder, f"fut_prof_{name}", options, settings.met_version)
    return folder, day_folders


def split_plan(plan, part, n_parts):
    """Release groups of part ``part`` (1..n_parts) of a plan cut into ``n_parts`` runs.

    The groups are taken in plan order (release time) and cut into consecutive chunks with about
    equal numbers of releases; a group is never divided. The part's releases are renumbered
    SourceID_1..n (NAME needs them in order); ``day_source_name`` keeps the name in the whole
    day's plan.

    Raises:
        ValueError if ``part`` is out of range or the part would have no release group.
    """
    if not 1 <= part <= n_parts:
        raise ValueError(f"Part {part} of {n_parts}?")
    groups = list(dict.fromkeys(plan["release_group"]))  # in plan order
    sizes = plan.groupby("release_group", sort=False).size().reindex(groups).to_numpy()
    # A group goes to the part its middle release falls in, counting releases in plan order.
    middle = np.cumsum(sizes) - sizes / 2
    part_of_group = 1 + np.floor(middle * n_parts / len(plan)).astype(int)  # 1..n_parts
    chosen = {g for g, k in zip(groups, part_of_group, strict=True) if k == part}
    if not chosen:
        raise ValueError(f"Part {part} of {n_parts} has no release group")
    subset = plan[plan["release_group"].isin(chosen)].copy()
    if "day_source_name" not in subset:  # a multi-day plan has it already (combine_release_plans)
        subset.insert(0, "day_source_name", subset["source_name"])
    subset["source_name"] = [f"SourceID_{n}" for n in range(1, len(subset) + 1)]
    return subset.reset_index(drop=True)


def prepare_split(month_folder, days, n_parts, parts=None, changes=None, label=None):
    """Runs ``profile_<label>_part<k>`` sharing the release groups of one or more prepared days.

    days: the days (consecutive days are combined in one plan first, as ``batch3days``).
    n_parts: how many runs the groups are cut into (``split_plan``); parts: which of them to
    prepare (default: all). changes: Options to change from the first day's (e.g. threads).
    label: the folder name (default ``split<N>_<first day>[-<last day>]``).
    Each run reuses the days' observation files and surface-pressure runs (already done);
    ``write_main_input`` keeps only its part. Returns the folders.
    """
    month_folder = Path(month_folder)
    days = [pd.Timestamp(d) for d in days]
    day_folders = [month_folder / f"{d:%Y%m%d}" for d in days]
    day_info = json.loads((day_folders[0] / "day.json").read_text())
    changes = changes or {}
    options = dataclasses.replace(Options(**day_info["options"]), **changes)
    if label is None:
        span = f"{days[0]:%Y%m%d}" + (f"-{days[-1]:%Y%m%d}" if len(days) > 1 else "")
        label = f"split{n_parts}_{span}"
    settings = run_settings(SITES[day_info["site"]], month_folder, days[0], days[-1], options)
    folders = []
    for part in parts or range(1, n_parts + 1):
        folder = month_folder / f"profile_{label}_part{part}"
        (folder / "output").mkdir(parents=True)
        info = {**day_info, "options": dataclasses.asdict(options),
                "day_folders": [str(f) for f in day_folders], "variant": f"{label}_part{part}",
                "changes": changes, "part": [part, n_parts]}  # fmt: skip
        (folder / "day.json").write_text(json.dumps(info, indent=2) + "\n")
        write_main_job(folder, f"fut_{label}_p{part}", options, settings.met_version)
        folders.append(folder)
        print(folder)
    return folders


def prepare_month(product, site_name, month, work, options, profiling=False, days=None):
    """Every day of ``month`` (or ``days``), the met-restore job, profiling variants, submit.sh."""
    site = SITES[site_name]
    month = pd.Period(month, freq="M")
    month_folder = Path(work) / f"{site.name}_{product.replace(' ', '_')}_{month:%Y%m}_future"
    if month_folder.exists():
        raise SystemExit(f"{month_folder} already exists; choose a new folder")
    month_folder.mkdir(parents=True)
    all_days = pd.date_range(month.start_time, month.end_time.normalize(), freq="D")
    if days is not None:
        all_days = pd.DatetimeIndex(pd.to_datetime(list(days))).normalize()
    summary = pd.DataFrame([prepare_day(product, site, d, month_folder, options)
                            for d in all_days])  # fmt: skip
    summary.to_csv(month_folder / "summary.csv", index=False)
    prepared = summary[summary["status"] == "prepared"]
    if prepared.empty:
        raise SystemExit("No day had soundings to run")
    versions = prepared["met_version"].unique()
    if len(versions) != 1:
        raise ValueError(f"The month spans met versions {list(versions)}: not supported yet")
    write_met_restore_job(month_folder / "restore_met.sbatch", month_folder, versions[0],
                          prepared["met_start"].min(), prepared["met_end"].max(),
                          MET_ROOT / versions[0], time="12:00:00")  # fmt: skip

    lines = ["#!/bin/bash", "# Met restore; per day the surface-pressure run, then the main run;",
             "# profiling runs after the surface-pressure runs of their days (afterok).",
             "set -e", f"cd {month_folder}", "restore=$(sbatch --parsable restore_met.sbatch)",
             'echo "met restore: job $restore"']  # fmt: skip
    for day in prepared["day"]:
        d = pd.Timestamp(day).strftime("%Y%m%d")
        lines += [f"ps_{d}=$(sbatch --parsable --dependency=afterok:$restore "
                  f"{d}/surface_pressure/job.sbatch)",
                  f"main=$(sbatch --parsable --dependency=afterok:$ps_{d} {d}/job.sbatch)",
                  f'echo "{d}: surface pressure $ps_{d}, main $main"']  # fmt: skip
    if profiling:
        for name, (variant_days, changes) in PROFILING_VARIANTS.items():
            if not all(day in set(prepared["day"]) for day in variant_days):
                print(f"Profiling variant {name} skipped: its days are not all prepared")
                continue
            folder, day_folders = prepare_variant(name, variant_days, changes, month_folder,
                                                  product, site, options)  # fmt: skip
            after = ":".join(f"$ps_{f.name}" for f in day_folders)
            lines += [f"job=$(sbatch --parsable --dependency=afterok:{after} "
                      f"{folder.name}/job.sbatch)", f'echo "profile_{name}: $job"']  # fmt: skip
    (month_folder / "submit.sh").write_text("\n".join(lines) + "\n")
    print(summary.to_string(index=False))
    print(f"Prepared {len(prepared)} days in {month_folder}; submit with: bash submit.sh")
    return summary


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if args and args[0] == "main-input":
        write_main_input(args[1])
    elif args and args[0] == "split":
        flags = dict(a.removeprefix("--").split("=", 1) for a in sys.argv[2:] if "=" in a)
        changes = {}
        if "threads" in flags:
            changes["threads"] = int(flags["threads"])
        if "output-grids" in flags:
            changes["output_grids"] = flags["output-grids"]
        parts = [int(k) for k in flags["parts"].split(",")] if "parts" in flags else None
        prepare_split(args[1], args[2].split(","), int(args[3]), parts, changes,
                      flags.get("label"))  # fmt: skip
    else:
        flags = sys.argv[1:]
        chosen = [a.removeprefix("--days=") for a in flags if a.startswith("--days=")]
        prepare_month(*args[:4], options=Options(test="--test" in flags),
                      profiling="--profiling" in flags,
                      days=chosen[0].split(",") if chosen else None)  # fmt: skip
