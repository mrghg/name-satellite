"""Small NAME tests for the future method: unbounded domain and nested output grids.

Part of the scoping work on the ``future-method`` branch (docs/future/scoping.md). Not part of the
workflow: nothing in src/ uses this.

The question: can NAME back-runs drop the side boundaries of the dispersion domain (particles run
for the full travel time, their 3D end positions give the boundary-condition term) and write
footprints on nested grids (fine near each release group, medium over the region, coarse
globally)? This script prepares a few small runs from an existing prepared day, so that run time,
memory and output size can be compared without running a whole day.

Releases: three release groups of the 1 Jan 2013 GOSAT-BRAZIL test day (16 levels each, 48 NAME
sources, ~73,000 particles, ~14% of the day), chosen far apart: NE coast, the Amazon near the
equator, and northern Argentina. They are copied unchanged from the day's ``release_plan.csv``.

Variants (one NAME run each, identical releases, patched NAME for the fate records, X20):

    bounded         today's domain and output grid (the baseline)
    shared_region   no side boundaries; one 0.1° grid over the whole region for every source,
                    plus the medium and coarse grids
    shared_global   as shared_region, but the 0.1° grid covers the globe (smoke test only: a
                    memory check of "memory grows with grid cells × sources")
    per_group       no side boundaries; a 0.1° grid of ±5° around each release group, used only
                    by that group's sources, plus the medium and coarse grids
    per_group_240h  as per_group, with 240 h maximum travel instead of 720 h

Grids (all aligned to whole degrees, so finer cells sum exactly into coarser ones):

    fine    0.1°  (per group ±5°, the region, or the globe)
    medium  0.5°  over today's SOUTHAMERICA dispersion domain, widened to whole degrees
    coarse  2.5°  globe

The top of the domain stays at 20,000 m agl (D21). Release met is left out (X21 is not under test).

Two stages:

    --smoke   every variant with 1 hour of travel, debug queue: checks that NAME accepts the
              input (an unbounded domain is not used anywhere else in our runs) and measures
              memory, which is set by the met and the grids from the start. Minutes.
    (default) bounded, per_group, per_group_240h with full travel, 4 threads each, run side by
              side: run time, particle-hours, output size and the end positions. Under an hour.

Usage (sci node):

    uv run python dev/future/prepare_unbounded_test.py DAY_FOLDER WORK [--smoke]

DAY_FOLDER is a prepared day folder with ``release_plan.csv`` (``release_group`` column), e.g.
``/work/scratch-pw5/mattrigby/name-satellite/rename_check_20261113/
GOSAT-BRAZIL_GOSAT-1_v10.0_201301/20130101``. WORK must not exist. Submit with
``bash WORK/submit.sh``; each variant is a subfolder with its own job.
"""

import dataclasses
import math
import sys
from pathlib import Path

import pandas as pd

from name_satellite import name_input
from name_satellite.domains import SITES
from name_satellite.lotus import write_job_script
from name_satellite.met import PATCHED_NAME_DIR
from name_satellite.prepare import run_settings

SITE = SITES["GOSAT-BRAZIL"]
DAY = pd.Timestamp("2013-01-01")
GROUPS = (
    "2013010114580250322002",  # NE coast of Brazil, 10.4°S
    "2013010116330280300018",  # Amazon, 0.2°N
    "2013010116420280360019",  # northern Argentina, 34.6°S
)
THREADS = 4
FINE_DEG, MEDIUM_DEG, COARSE_DEG = 0.1, 0.5, 2.5
FINE_HALF_WIDTH_DEG = 5.0  # per-group fine grid: ±5° around the group's release
# Medium grid and shared regional fine grid: today's dispersion domain, widened to whole degrees.
REGION = {"x_min": -97.0, "x_max": -19.0, "y_min": -65.0, "y_max": 28.0}
GLOBE = {"x_min": -180.0, "x_max": 180.0, "y_min": -90.0, "y_max": 90.0}


@dataclasses.dataclass(frozen=True)
class Grid:
    """A named Lat-Long output grid given by its outer edges and cell sizes (degrees)."""

    name: str
    dx: float
    dy: float
    x_min: float
    x_max: float
    y_min: float
    y_max: float

    @property
    def nx(self):
        return round((self.x_max - self.x_min) / self.dx)

    @property
    def ny(self):
        return round((self.y_max - self.y_min) / self.dy)

    def input_line(self):
        # NAME's X Min / Y Min are the centre of the first cell.
        centre_x = self.x_min + self.dx / 2
        centre_y = self.y_min + self.dy / 2
        return (
            f"{self.name}, Lat-Long, {self.nx}, {self.ny}, {self.dx}, {self.dy}, "
            f"{centre_x:.4f}, {centre_y:.4f}"
        )


def box_around(name, lon, lat, half_width, dx):
    """Grid of cell size dx covering lon/lat ± half_width, with edges on multiples of dx."""
    x_min = math.floor((lon - half_width) / dx) * dx
    y_min = math.floor((lat - half_width) / dx) * dx
    width = 2 * half_width
    return Grid(name, dx, dx, round(x_min, 6), round(x_min + width, 6), round(y_min, 6),
                round(y_min + width, 6))  # fmt: skip


def variant_grids(variant, plan):
    """{source_name: [grids it is written on]} and the list of all grids, for one variant."""
    medium = Grid("HGridMedium", MEDIUM_DEG, MEDIUM_DEG, **REGION)
    coarse = Grid("HGridCoarse", COARSE_DEG, COARSE_DEG, **GLOBE)
    if variant == "bounded":
        g = SITE.domain.grid  # today's grid (cell centres in the site definition)
        today = Grid("HGrid1", g.dx, g.dy, g.x_min - g.dx / 2, g.x_min + (g.nx - 0.5) * g.dx,
                     g.y_min - g.dy / 2, g.y_min + (g.ny - 0.5) * g.dy)  # fmt: skip
        return {s: [today] for s in plan["source_name"]}, [today]
    if variant == "shared_region":
        fine = Grid("HGridFine", FINE_DEG, FINE_DEG, **REGION)
        grids = [fine, medium, coarse]
        return {s: grids for s in plan["source_name"]}, grids
    if variant == "shared_global":
        fine = Grid("HGridFine", FINE_DEG, FINE_DEG, **GLOBE)
        grids = [fine, medium, coarse]
        return {s: grids for s in plan["source_name"]}, grids
    if variant.startswith("per_group"):
        by_source, fines = {}, []
        for number, (_group, rows) in enumerate(plan.groupby("release_group", sort=True), 1):
            fine = box_around(f"HGridFineG{number}", rows["lon"].iloc[0], rows["lat"].iloc[0],
                              FINE_HALF_WIDTH_DEG, FINE_DEG)  # fmt: skip
            fines.append(fine)
            for source in rows["source_name"]:
                by_source[source] = [fine, medium, coarse]
        return by_source, [*fines, medium, coarse]
    raise ValueError(f"Unknown variant {variant!r}")


def grids_block(grids, settings, timing):
    """Horizontal, vertical and temporal grids (vertical and temporal as in name_input)."""
    lines = ["Horizontal Grids:", "Name, H-Coord, nX, nY, dX, dY, X Min, Y Min"]
    lines += [grid.input_line() for grid in grids]
    # Reuse name_input's vertical and temporal grids: everything after its horizontal grid.
    standard = name_input.grids_block(settings, timing)
    rest = standard[standard.index("Vertical Grids:") :]
    return "\n".join(lines) + "\n\n" + rest


def domain_block(variant, settings, timing):
    """Today's bounded domain, or no side boundaries (top kept at 20,000 m agl, D21).

    NAME refuses an H-Coord (or any horizontal bounds) for an unbounded domain.
    """
    if variant == "bounded":
        return name_input.domain_block(settings, timing)
    return (
        "Domains:\n"
        "Name, H Unbounded?, H-Coord, X Min, X Max, Y Min, Y Max, Z Unbounded?, Z-Coord, Z Max, "
        "T Unbounded?, Start Time, End Time, Max Travel Time\n"
        f"Dispersion Domain, Yes, , , , , , No, m agl, {settings.domain.z_max_m_agl}, "
        f"Yes, , , {timing['max_travel_hours']}:00\n"
    )


def sources_and_requests_block(plan, by_source, settings, timing):
    """Sources as in name_input; one field request per source and grid, one output group per grid.

    Each output group is one fields file (``<group>_C1_T1_<time>.txt``),
    e.g. ``Fields_HGridFineG2``.
    """
    standard = name_input.sources_and_requests_block(plan, settings, timing)
    sources = standard[: standard.index("Output Requirements - Fields:")]
    n_steps = timing["window_hours"] * 60 // settings.sync_minutes
    lines = [
        "Output Requirements - Fields:",
        "Name, Quantity, Species, Source, H-Grid, Z-Grid, T-Grid, BL Average, T Av Or Int, "
        "Av Time, # Av Times, Sync?, Output Route, Across, Separate File, Output Format, "
        "Output Group",
    ]
    for source in plan["source_name"]:
        number = source.removeprefix("SourceID_")
        for grid in by_source[source]:
            lines.append(
                f"Req_{grid.name}_{number}, Mixing Ratio, INERT-TRACER, {source}, {grid.name}, "
                f"ZGrid1, TGridEndOfRun, No, Int, {timing['window_hours']}:00, {n_steps}, No, D, "
                f"TZ, T, IA, Fields_{grid.name}"
            )
    return sources + "\n".join(lines) + "\n"


def input_text(variant, plan, settings):
    """The whole NAME input file for one variant (blocks as in name_input.name_input_text)."""
    timing = name_input.run_timing(plan, settings)
    by_source, grids = variant_grids(variant, plan)
    blocks = [
        f"! NAME input file: future-method test, variant {variant} "
        "(dev/future/prepare_unbounded_test.py)\n",
        name_input.main_options_block(settings, len(plan)),
        name_input.restart_and_cases_block(),
        name_input.openmp_block(settings),
        name_input.timer_block(settings),
        name_input.files_and_coordinates_block(settings),
        grids_block(grids, settings, timing),
        domain_block(variant, settings, timing),
        name_input.species_block(),
        name_input.dispersion_options_block(settings),
        name_input.particle_fate_block(settings),
        sources_and_requests_block(plan, by_source, settings, timing),
        name_input.met_declarations_block(settings),
    ]
    return "\n".join(blocks)


def add_memory_sampling(job_path, input_path):
    """Make the job record NAME's peak memory (kB) in ``peak_memory_kb.txt``.

    Slurm samples memory too rarely for the short smoke runs, and timing ``mpirun`` measures
    mpirun, not NAME. So NAME runs in the background while the job reads the peak resident memory
    (VmHWM, kept by Linux for each process) of the processes whose command line contains this
    run's input file, every 2 s. Several variants may share a node; the input path tells them apart.
    """
    job = job_path.read_text()
    run_line = next(line for line in job.splitlines() if line.startswith("time mpirun"))
    sampled = f"""{run_line.removeprefix("time ")} &
name_pid=$!
peak=0
while kill -0 $name_pid 2>/dev/null; do
    for p in $(pgrep -f "{input_path}"); do
        hwm=$(awk '/VmHWM/ {{print $2}}' /proc/$p/status 2>/dev/null)
        if [ -n "$hwm" ] && [ "$hwm" -gt "$peak" ]; then peak=$hwm; fi
    done
    echo $peak > peak_memory_kb.txt
    sleep 2
done
wait $name_pid"""
    if job.count(run_line) != 1 or "status=$?" not in job:
        raise ValueError(f"Unexpected job script layout in {job_path}")
    job_path.write_text(job.replace(run_line, sampled))


def main(day_folder, work, smoke):
    day_folder, work = Path(day_folder), Path(work)
    if work.exists():
        raise SystemExit(f"{work} already exists; choose a new folder")
    plan = pd.read_csv(day_folder / "release_plan.csv", dtype={"release_group": str},
                       parse_dates=["release_start", "release_end"])  # fmt: skip
    plan = plan[plan["release_group"].isin(GROUPS)].copy()
    if plan["release_group"].nunique() != len(GROUPS):
        raise ValueError(f"Not all of {GROUPS} are in {day_folder}/release_plan.csv")
    # Renumber the sources 1..n (NAME's source numbers in the fate records follow input order).
    plan["source_name"] = [f"SourceID_{n}" for n in range(1, len(plan) + 1)]

    if smoke:
        variants = {name: 1 for name in ("bounded", "shared_region", "shared_global", "per_group")}
        qos, time = "debug", "00:30:00"
    else:
        variants = {"bounded": 720, "per_group": 720, "per_group_240h": 240}
        qos, time = "high", "04:00:00"

    submit = ["#!/bin/bash", "# Submit every variant; they are independent.", "set -e"]
    for variant, hours in variants.items():
        folder = work / variant
        (folder / "output").mkdir(parents=True)
        settings = run_settings(SITE, folder, DAY, THREADS, test=False)
        settings = dataclasses.replace(settings, run_name=variant, max_age_hours=hours,
                                       release_met=False)  # fmt: skip
        (folder / "input.txt").write_text(input_text(variant, plan, settings))
        plan.to_csv(folder / "release_plan.csv", index=False)
        write_job_script(folder / "job.sbatch", folder, folder / "input.txt", threads=THREADS,
                         qos=qos, time=time, memory="32G", job_name=f"fut_{variant}",
                         name_dir=PATCHED_NAME_DIR)  # fmt: skip
        add_memory_sampling(folder / "job.sbatch", folder / "input.txt")
        submit.append(f"sbatch --parsable {folder / 'job.sbatch'}")
    (work / "submit.sh").write_text("\n".join(submit) + "\n")
    particles = int((plan["particles_per_second"] * 60 + 1).sum())
    print(f"Prepared {work}: {', '.join(variants)}; {len(plan)} sources, {particles} particles "
          f"per run. Submit with: bash {work / 'submit.sh'}")  # fmt: skip


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], smoke="--smoke" in sys.argv[3:])
