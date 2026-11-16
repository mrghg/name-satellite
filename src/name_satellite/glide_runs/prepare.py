"""Write one GLIDE reference run: release plan, NAME input file and Slurm job.

Where this sits
---------------
    spec.py, plan.py  ->  **prepare.py: run folder**  ->  LOTUS2 (NAME)  ->  deliver.py (§8 files)

A run is one receptor, a list of hourly release times, and a variant of the settings (``main``,
``seed1`` ... ``seed5``, ``nomeso``, ``group``; spec.VARIANTS). It is written to its own folder:

* ``plan.csv``: the release plan (``plan.point_release_plan``);
* ``input.txt``: the NAME input file. Besides the footprints and particle-fate records, it asks
  for met at the receptor (``ReceptorMet_<R...>_C1.txt``: model ground height, land fraction,
  surface pressure and boundary-layer depth at the first release time), for the run record
  (§3, §8). With ``positions_hours``, it also asks for every particle's position on an hourly
  grid (``Positions_C1.txt``; single-release runs, §9);
* ``job.sbatch``: runs NAME (the ACRG-patched build by default, whose fate records name each
  particle and why it stopped) and fails if NAME's error file is not empty;
* ``run.json``: what was asked for (receptor, times, variant, settings) for the run record.

The met (Mk10, plain files: restoring makes links) is restored by a separate job,
``restore_met.sbatch``, into ``met_folder``, which several runs may share; submit it first and
make NAME jobs depend on it.

Command line (from the repository folder)::

    uv run python -m name_satellite.glide_runs.prepare RUN_FOLDER RECEPTOR FIRST LAST \\
        [--variant main] [--max-age 240] [--positions-hours N] [--threads 4] [--qos high] \\
        [--time 12:00:00] [--memory 32G] [--met-folder DIR] [--unpatched] \\
        [--surface-pressure CSV]   # Lamont only: from glide_runs.column

FIRST and LAST are the first and last release start times (UTC, e.g. 2019-01-08T00).
``--max-age`` shortens particle travel for tests only (the spec fixes 240 h).
"""

import argparse
import json
from dataclasses import asdict, replace
from pathlib import Path

import pandas as pd

from name_satellite import lotus
from name_satellite.glide_runs import spec
from name_satellite.glide_runs.column import COLUMN_METHOD, column_release_plan
from name_satellite.glide_runs.plan import point_release_plan
from name_satellite.met import NAME_DIR, PATCHED_NAME_DIR
from name_satellite.name_input import (
    ParticlePositionOutput,
    name_input_text,
    run_timing,
)
from name_satellite.release_met import release_met_blocks

# Met requested at each receptor (release_met: suffix, NAME quantity, level). "ground" is 0 m
# above the model's ground. Topography is the model's ground height above sea level (m), Land
# Fraction the model's land fraction (0-1) of that point (NAME documentation A01, output
# quantities); both answer §3's "what the model sees at each receptor".
RECEPTOR_MET = [
    ("TOPO", "Topography", "ground"),
    ("LAND", "Land Fraction", "ground"),
    ("PS", "Pressure (Pa)", "ground"),
    ("BL", "Boundary Layer Depth", "ground"),
]


def receptor_met_text(plan, grid):
    """NAME input requesting the met of ``RECEPTOR_MET`` at the run's receptor.

    Sampled once, at the first release time. Output group ``ReleaseMet_R1``
    (file ``ReleaseMet_R1_C1.txt``; read with ``release_met.read_release_met_file``).
    """
    first = plan.iloc[0]
    points = pd.DataFrame(
        {
            "sounding": ["R1"],
            "lon": [first["lon"]],
            "lat": [first["lat"]],
            "time": [first["release_time"]],
            "group": [1],
        }
    )
    return release_met_blocks(points, grid.dx, grid.dy, quantities=RECEPTOR_MET)


def run_settings(run_folder, receptor, variant, max_age_hours, positions, threads, met_folder):
    """``RunSettings`` for this run: spec.reference_settings with the variant's changes."""
    domain = spec.DOMAINS[receptor.domain]
    max_particles = spec.MAX_PARTICLES_COLUMN if receptor.is_column else spec.MAX_PARTICLES_POINT
    settings = spec.reference_settings(
        run_name=Path(run_folder).name,
        output_folder=Path(run_folder) / "output",
        domain=domain,
        met_folder=met_folder,
        max_particles=max_particles,
        threads=threads,
    )
    settings = spec.variant_settings(settings, variant)
    return replace(settings, max_age_hours=max_age_hours, particle_positions=positions)


def write_run_folder(
    run_folder, plan, settings, job_name, qos, time, memory, patched=True, record=None
):
    """Write ``plan.csv``, ``input.txt``, ``job.sbatch``, ``restore_met.sbatch`` and ``run.json``.

    Args:
        run_folder: the run's folder (created here; it must not exist).
        plan: release plan (``plan.py``; any plan ``name_input`` accepts).
        settings: ``RunSettings`` (threads, met folder and output folder are taken from it).
        job_name, qos, time, memory: Slurm job (``lotus.write_job_script``).
        patched: run the ACRG-patched NAME (fate records with particle index and stop reason).
        record: extra entries for ``run.json``.

    The input file also requests the receptor met (``receptor_met_text``) at the plan's first
    release. The met restore covers the integration window plus an hour either side.
    """
    run_folder = Path(run_folder)
    run_folder.mkdir(parents=True)
    (run_folder / "output").mkdir()
    plan.to_csv(run_folder / "plan.csv", index=False)
    text = name_input_text(plan, settings) + "\n" + receptor_met_text(plan, settings.grid)
    (run_folder / "input.txt").write_text(text)

    timing = run_timing(plan, settings)
    lotus.write_job_script(
        run_folder / "job.sbatch",
        run_folder,
        run_folder / "input.txt",
        threads=settings.threads,
        qos=qos,
        time=time,
        memory=memory,
        job_name=job_name,
        name_dir=PATCHED_NAME_DIR if patched else NAME_DIR,
    )
    # Met for the whole integration window, plus an hour either side for interpolation.
    met_start = timing["window_start"] - pd.Timedelta(hours=1)
    met_end = plan["release_start"].max() + pd.Timedelta(hours=1)
    lotus.write_met_restore_job(
        run_folder / "restore_met.sbatch",
        run_folder,
        settings.met_version.name,
        met_start.strftime("%Y-%m-%dT%H:%M"),
        met_end.strftime("%Y-%m-%dT%H:%M"),
        settings.met_folder,
        time="02:00:00",
    )

    full_record = {
        "spec_version": spec.SPEC_VERSION,
        **(record or {}),
        "n_releases": len(plan),
        "name_build": "ACRG-patched 8.7.1 (fate records)" if patched else "NAME 8.7.1",
        "met_period": [str(met_start), str(met_end)],
        "settings": {
            key: str(value) if not isinstance(value, (int, float, bool, type(None))) else value
            for key, value in asdict(settings).items()
        },
    }
    (run_folder / "run.json").write_text(json.dumps(full_record, indent=2, default=str) + "\n")
    return run_folder


def prepare_run(
    run_folder,
    receptor_code,
    first_release,
    last_release,
    variant="main",
    max_age_hours=spec.MAX_AGE_HOURS,
    positions_hours=None,
    threads=4,
    qos="high",
    time="12:00:00",
    memory="32G",
    met_folder=None,
    patched=True,
    surface_pressure_csv=None,
):
    """Write a run folder (see the module docstring). Returns the folder.

    Args:
        run_folder: new folder for the run (must not exist: one folder per attempt).
        receptor_code: a key of ``spec.RECEPTORS``.
        first_release, last_release: first and last release start times (UTC, hourly).
        variant: a key of ``spec.VARIANTS``.
        max_age_hours: particle travel time (h); only tests change it from 240.
        positions_hours: if given, also write every particle's position each hour, for this
            many hours back from the end of the last release (single-release runs: 241).
        threads, qos, time, memory: Slurm job (``lotus.write_job_script``).
        met_folder: where the met is restored (default ``<run_folder>/met``).
        patched: run the ACRG-patched NAME (fate records with particle index and stop reason).
        surface_pressure_csv: for the column receptor (Lamont) only: the model surface pressure
            at each release start (``column.read_surface_pressure_runs``, saved as CSV); it must
            cover every release of the run.

    Raises:
        FileExistsError if ``run_folder`` exists; ValueError for unknown receptors or variants,
        or a column run without (complete) surface pressures.
    """
    run_folder = Path(run_folder)
    if run_folder.exists():
        raise FileExistsError(f"{run_folder} exists: use a new folder for each attempt")
    receptor = spec.RECEPTORS[receptor_code]
    times = pd.date_range(first_release, last_release, freq="1h")
    if receptor.is_column:
        if surface_pressure_csv is None:
            raise ValueError(f"{receptor_code} is a column: give its surface pressures (§5)")
        pressures = pd.read_csv(surface_pressure_csv, parse_dates=["release_time"])
        pressures = pressures[pressures["release_time"].isin(times)]
        if len(pressures) != len(times):
            raise ValueError(f"{surface_pressure_csv} lacks some of the run's release times")
        plan = column_release_plan(receptor, pressures)
    else:
        plan = point_release_plan(receptor, times)
    met_folder = Path(met_folder) if met_folder else run_folder / "met"

    positions = None
    if positions_hours is not None:
        # T0 = end of the last release (the chronologically latest time); NAME counts a
        # backwards run's time grid back from T0.
        positions = ParticlePositionOutput(
            first_time=plan["release_start"].max(), n_times=positions_hours, interval_minutes=60
        )
    settings = run_settings(
        run_folder, receptor, variant, max_age_hours, positions, threads, met_folder
    )

    record = {
        "receptor": asdict(receptor),
        "first_release": str(times[0]),
        "last_release": str(times[-1]),
        "variant": variant,
    }
    if receptor.is_column:
        record["column_method"] = COLUMN_METHOD
    write_run_folder(
        run_folder,
        plan,
        settings,
        job_name=f"glide_{receptor_code}_{variant}",
        qos=qos,
        time=time,
        memory=memory,
        patched=patched,
        record=record,
    )
    return run_folder


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("run_folder")
    parser.add_argument("receptor", choices=sorted(spec.RECEPTORS))
    parser.add_argument("first_release")
    parser.add_argument("last_release")
    parser.add_argument("--variant", default="main", choices=sorted(spec.VARIANTS))
    parser.add_argument("--max-age", type=int, default=spec.MAX_AGE_HOURS)
    parser.add_argument("--positions-hours", type=int, default=None)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--qos", default="high")
    parser.add_argument("--time", default="12:00:00")
    parser.add_argument("--memory", default="32G")
    parser.add_argument("--met-folder", default=None)
    parser.add_argument("--unpatched", action="store_true")
    parser.add_argument("--surface-pressure", default=None, help="CSV (column receptor only)")
    args = parser.parse_args(argv)
    folder = prepare_run(
        args.run_folder,
        args.receptor,
        args.first_release,
        args.last_release,
        variant=args.variant,
        max_age_hours=args.max_age,
        positions_hours=args.positions_hours,
        threads=args.threads,
        qos=args.qos,
        time=args.time,
        memory=args.memory,
        met_folder=args.met_folder,
        patched=not args.unpatched,
        surface_pressure_csv=args.surface_pressure,
    )
    print(f"Wrote {folder}")


if __name__ == "__main__":
    main()
