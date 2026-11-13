"""Prepare a month of NAME runs from Level-2 satellite data: workflow steps 1-4 chained together.

Where this sits (see docs/workflow.md)
--------------------------------------
    L2 files (CEDA)  ->  **prepare**: per day, observation file + surface-pressure run + job
    ->  on LOTUS2: met restore -> surface-pressure run -> main job (writes the NAME input, runs
    NAME)  ->  read and process the output (steps 5-6)

Each day has TWO NAME runs, because the release layers are scaled to the model's surface pressure
(decision D19), which only a NAME run can report:

``prepare_day`` (now, on the sci node):
  1. downloads the day's L2 file (``l2_files``; days without a file are skipped and reported);
  2. reads it and keeps good soundings in the site's box (``observations``: quality flag, defined
     pressure levels, pressures decreasing from the surface, lat/lon bounds); legacy binning only
     if the site asks for it (validation only; none for GOSAT-BRAZIL);
  3. writes the observation file, keyed by sounding ID (D18), and ``day.json`` (what was asked);
  4. writes the met-only surface-pressure run (``surface_pressure``: model surface pressure and
     domain-top pressure at each sounding) and its job;
  5. writes the main job, which first runs ``write_main_input`` and then NAME.

``write_main_input`` (in the main job, once the surface-pressure run has finished):
  6. reads the model surface pressure and checks it against the retrieval's (model / retrieval
     within 0.9-1.1; ``surface_pressure.csv``);
  7. scales every sounding's levels to its model surface pressure (D19); checks that the
     modelled levels (``releases.MODELLED_LEVELS``, 16 for GOSAT) lie inside the domain and carry
     at least 0.85 of the column weight (D21); groups the soundings that share releases
     (``release_groups``, D20; ``release_groups.csv``) and plans one set of releases per group
     (``releases.build_release_plan``, ``release_plan.csv``), with each level's particle rate
     from ``releases.PARTICLE_RATES`` (D22);
  8. writes the NAME input ``input.txt`` (with release-point met, X21). The job then runs NAME.

``prepare_month`` also writes one met-restore job for the whole month and ``submit.sh``, which
submits the restore job, then each day's surface-pressure job after it and the day's main job
after that (Slurm ``afterok``: a job starts only if the one before succeeded).

Command line
------------
::

    uv run python -m name_satellite.prepare PRODUCT SITE YYYY-MM WORK [--test] [--days=D1,D2,...]
    uv run python -m name_satellite.prepare "GOSAT-1 v10.0" GOSAT-BRAZIL 2013-01 WORK
    uv run python -m name_satellite.prepare main-input WORK/<month folder>/YYYYMMDD

PRODUCT is a key of ``l2_files.PRODUCTS``, SITE of ``domains.SITES``. The month folder must not
exist yet. ``--test``: particles travel 1 hour only, debug queue, to check a set-up quickly
(lesson of X20); ``processing`` refuses such runs. ``--days=2013-01-01,2013-01-02``: prepare only
these days of the month. ``main-input`` is what each day's main job runs; it can be run by hand
once the surface-pressure run has finished (it overwrites the files it writes). Then submit with
``bash submit.sh`` in the month folder.

Files (month folder ``WORK/<SITE>_<PRODUCT, spaces as _>_<YYYYMM>/``, e.g.
``GOSAT-BRAZIL_GOSAT-1_v10.0_201301/``)
-------------------------------------------------------------------------------------------------
Written now, by ``prepare_month``:

* ``summary.csv``: one row per day asked for: ``day``, ``status`` ("prepared", "no L2 file",
  "no soundings in the box"), ``l2_soundings`` (in the L2 file), ``soundings`` (selected),
  ``release_groups``, ``sources`` (NAME releases), ``met_version``, ``met_start``/``met_end``
  (period of met the day needs). Days without soundings get no day folder.
* ``restore_met.sbatch``: Slurm job restoring the month's met (``met_restore``); logs
  ``restore-<jobid>.out/.err``.
* ``submit.sh``: submits the restore job, then each day's surface-pressure job after it, then the
  day's main job after that (Slurm ``afterok``: a job starts only if the one before succeeded).
* ``YYYYMMDD/obs_<SITE>_<YYYYMMDD>.nc``: the observation file (``observations.to_obs``; netCDF, one
  entry per sounding, dim ``sounding_id``). Read by ``write_main_input`` and ``processing``.
* ``YYYYMMDD/day.json``: product, site, day, threads, test: what was asked. Read by
  ``write_main_input`` and ``processing``.
* ``YYYYMMDD/job.sbatch``: the main job (patched NAME, ``met.PATCHED_NAME_DIR``, X20).
* ``YYYYMMDD/surface_pressure/soundings.csv``: one row per sounding: ``sounding`` (label S0001,
  ...), ``point_number``, ``sounding_id``, ``lon``, ``lat`` (degrees), ``time`` (UTC, whole
  seconds), ``z_pa`` (centre of its unscaled level-1 layer, Pa; not used by this run), ``group``
  (output group of the location requests). Read back by ``write_main_input``.
* ``YYYYMMDD/surface_pressure/input.txt`` and ``job.sbatch``: the met-only NAME run (stock NAME).

Written on LOTUS2:

* ``YYYYMMDD/surface_pressure/output/ReleaseMet_S0001_C1.txt``, ...: one NAME point-output file
  per sounding (surface pressure and domain-top pressure, Pa; ``surface_pressure``).
* By ``write_main_input`` (in the main job), in ``YYYYMMDD/``:

  - ``surface_pressure.csv``: one row per sounding: ``sounding_id``, ``retrieval_hpa`` (level 0),
    ``model_hpa``, ``ratio`` (model / retrieval, applied to every level, D19), ``domain_top_hpa``
    (model pressure at the 20 km agl domain top, D21);
  - ``release_groups.csv``: one row per sounding: ``sounding_id``, ``release_group`` (the group's
    ID = sounding ID of its earliest member, D20);
  - ``release_plan.csv``: one row per NAME release (source), columns as
    ``releases.build_release_plan`` (``source_name`` SourceID_<n>, ``point_number`` = group
    number, ``level_number`` 1-16, ``release_end``/``release_start`` UTC, ``lon``/``lat``/
    ``dlon``/``dlat`` degrees, ``z_pa``/``dz_pa`` Pa, ``particles_per_second``, and
    ``release_group`` = the group's ID, i.e. its earliest sounding's ID);
  - ``input.txt``: the main NAME input (``name_input``).

* ``YYYYMMDD/output/``: the main run's NAME output (``Fields_grid1*_C1_T1_*.txt``,
  ``Particles_C1.txt``, ``ReleaseMet_S0001_C1.txt``, ...; see ``name_output``), read by
  ``processing``.
* Next to each ``input.txt``: NAME's ``inputLog.txt``, ``inputError.txt`` (must be empty; the job
  fails otherwise) and ``inputTiming.txt``; and Slurm logs ``slurm-<jobid>.out/.err``.

To check a day: ``head`` the CSV files; ``ncdump -h`` the observation file; ``grep -c SourceID
input.txt`` counts sources (twice: one source and one field request per release).

NAME versions: the main run uses the ACRG-patched NAME 8.7.1 (``met.PATCHED_NAME_DIR``), whose
fate records also carry each particle's index and fate code (X20). The fate records have no met
columns (heights are m above ground, D24, so no pressure is needed). The met-only run uses the
stock NAME 8.7.1.
"""

import dataclasses
import getpass
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from name_satellite.domains import SITES
from name_satellite.l2_files import L2FileMissing, download_l2
from name_satellite.lotus import REPO, write_job_script, write_met_restore_job
from name_satellite.met import PATCHED_NAME_DIR, met_version_for_period
from name_satellite.name_input import RunSettings, run_timing, write_name_input
from name_satellite.observations import (
    bad_pressure_mask,
    latlon_mask,
    legacy_bin,
    level_order_mask,
    quality_mask,
    read_l2,
    select,
    to_obs,
    write_obs_file,
)
from name_satellite.release_groups import group_soundings, release_points
from name_satellite.release_met import MAX_SOUNDINGS, release_met_soundings
from name_satellite.releases import (
    MODELLED_LEVELS,
    PARTICLE_RATES,
    build_release_plan,
    check_modelled_levels,
    scale_to_model_surface,
)
from name_satellite.surface_pressure import (
    read_surface_pressure,
    surface_pressure_input_text,
    surface_pressure_settings,
    surface_pressure_table,
)

# Working folders in scratch, per user. Set NAME_SATELLITE_SCRATCH to move them, and
# NAME_SATELLITE_MET to use met restored elsewhere (e.g. someone else's readable copy for the same
# period: about 550 GB for a month of Mk6). Slurm jobs inherit these from the submitting shell.
SCRATCH = Path(
    os.environ.get(
        "NAME_SATELLITE_SCRATCH", f"/work/scratch-pw5/{getpass.getuser()}/name-satellite"
    )
)
L2_CACHE = SCRATCH / "l2_cache"
MET_ROOT = Path(os.environ.get("NAME_SATELLITE_MET", SCRATCH / "met"))
THREADS = 4  # about 80% efficient (benchmark, dev/notes/name-behaviour.md)


def select_soundings(table, site):
    """The good soundings of one day's table in the site's box (binned if the site says so).

    Good = quality flag 0, all 20 pressure levels defined, pressure decreasing from the surface
    (``observations``). Box: lower bounds included, upper excluded (degrees).

    Args:
        table: ``observations.read_l2`` output (dims sounding, level).
        site: ``domains.Site``.

    Returns:
        The selected sounding table (same layout; legacy-binned if ``site.bin_size`` is set).
    """
    good = quality_mask(table) & bad_pressure_mask(table) & level_order_mask(table)
    table = select(table, good)
    if site.bin_size is not None:
        # legacy_bin applies the box itself (soundings outside the bins are dropped).
        return legacy_bin(
            table, site.lat_bounds, site.lon_bounds, site.bin_size, footprint_km=site.footprint_km
        )
    return select(table, latlon_mask(table, site.lat_bounds, site.lon_bounds))


def run_settings(site, run_folder, day, threads, test):
    """NAME settings for one day's run (legacy defaults apart from the deviations adopted).

    The met version must cover from 31 days before ``day`` to the day after (30-day travel plus
    margin); output goes to ``run_folder/output``, met to ``MET_ROOT/<met version>``. With
    ``test`` the particles travel 1 hour only. ``processing.read_day`` calls this again to
    record the settings in the footprint file, so it must stay deterministic.

    Raises:
        ValueError (``met.met_version_for_period``) if that period spans a met-version change.
    """
    met_version = met_version_for_period(day - pd.Timedelta(days=31), day + pd.Timedelta(days=1))
    settings = RunSettings(
        run_name=run_folder.name,
        output_folder=str(run_folder / "output"),
        grid=site.domain.grid,
        domain=site.domain.dispersion,
        met_version=met_version,
        met_folder=str(MET_ROOT / met_version.name),
        threads=threads,
    )
    if test:
        settings = dataclasses.replace(settings, max_age_hours=1)
    return settings


def release_boxes(obs):
    """(dlon, dlat) per sounding for binned soundings, else None (the footprint box is used)."""
    return (obs["dlon"].values, obs["dlat"].values) if "dlat" in obs else None


def grouped_release_plan(obs, site, n_levels, pressure_levels_hpa=None, particle_rates=None):
    """Release groups (D20) and their release plan for one day.

    Args:
        obs: the day's observation Dataset (one entry per sounding).
        site: the ``domains.Site`` (release cells, footprint size, output grid).
        n_levels: number of modelled levels (``releases.MODELLED_LEVELS``, D21).
        pressure_levels_hpa: optional (sounding, level) levels to release at (D19-scaled).
        particle_rates: particles/s per level (``releases.PARTICLE_RATES``, D22); default legacy.

    Returns:
        (plan, members): the release plan (one set of releases per group, with a
        ``release_group`` column) and the sounding -> group table.
    """
    groups = group_soundings(obs, site.release_cell, site.domain.grid, site.release_group_seconds)
    points, members = release_points(obs, groups, site.footprint_km, pressure_levels_hpa)
    boxes = (points["dlon"].values, points["dlat"].values)
    plan = build_release_plan(points, n_levels=n_levels, footprint_diameter_km=site.footprint_km,
                              box_degrees=boxes, particle_rates=particle_rates)  # fmt: skip
    return plan, members


def prepare_day(product, site, day, month_folder, threads=THREADS, test=False):
    """Prepare one day: observation file, surface-pressure run and the main job.

    Writes the day folder ``month_folder/YYYYMMDD`` (see the module docstring), but only if the
    day has an L2 file and soundings in the box. The release plan made here uses the retrieval's
    own levels; only the run timing and the number of sources are taken from it (scaling the
    levels, D19, changes only the pressures). ``write_main_input`` makes the final plan.

    Args:
        product: L2 product name (``l2_files.PRODUCTS``), e.g. "GOSAT-1 v10.0".
        site: ``domains.Site``.
        day: the day (anything ``pandas.Timestamp`` reads).
        month_folder: the month's run folder (exists).
        threads: OpenMP threads (= CPUs) of the main run.
        test: 1-hour travel and the debug queue.

    Returns:
        A summary dict (one row of summary.csv).

    Raises:
        ValueError if the day has more than ``release_met.MAX_SOUNDINGS`` soundings (one met-only
        run cannot sample more), or a met-version change falls in the day's met period.
    """
    day = pd.Timestamp(day)
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
        # Release-point met allows at most MAX_SOUNDINGS per run (release_met.py).
        raise ValueError(f"{day:%Y-%m-%d}: {row['soundings']} soundings > {MAX_SOUNDINGS}")

    run_folder = month_folder / f"{day:%Y%m%d}"
    (run_folder / "output").mkdir(parents=True)
    write_obs_file(table, run_folder / f"obs_{site.name}_{day:%Y%m%d}.nc")
    info = {"product": product, "site": site.name, "day": f"{day:%Y-%m-%d}", "threads": threads,
            "test": test}  # fmt: skip
    (run_folder / "day.json").write_text(json.dumps(info, indent=2) + "\n")

    # Provisional plans from the retrieval's own levels: times and positions are final (only the
    # pressures change, D19). Per sounding: where the met-only run samples. Grouped (D20): the
    # number of releases and the run timing.
    obs = to_obs(table)
    n_levels = MODELLED_LEVELS[product]
    sounding_plan = build_release_plan(obs, n_levels=n_levels,
                                       footprint_diameter_km=site.footprint_km,
                                       box_degrees=release_boxes(obs))  # fmt: skip
    plan, members = grouped_release_plan(obs, site, n_levels,
                                         particle_rates=PARTICLE_RATES[product])  # fmt: skip
    settings = run_settings(site, run_folder, day, threads, test)
    qos = "debug" if test else "standard"

    # The met-only run: model surface pressure at each sounding.
    folder = run_folder / "surface_pressure"
    (folder / "output").mkdir(parents=True)
    soundings = release_met_soundings(sounding_plan)
    soundings.to_csv(folder / "soundings.csv", index=False)
    sp_settings = surface_pressure_settings(settings, folder / "output")
    (folder / "input.txt").write_text(surface_pressure_input_text(soundings, sp_settings))
    write_job_script(folder / "job.sbatch", folder, folder / "input.txt", threads=1, qos=qos,
                     time="00:30:00", memory="8G",
                     job_name=f"{site.name}_{day:%Y%m%d}_ps")  # fmt: skip

    # The main job: write the NAME input from the surface pressure, then run NAME.
    write_main = (
        f"cd {REPO}\n"
        f"~/.local/bin/uv run --no-sync python -m name_satellite.prepare main-input {run_folder}"
        " || exit 1"
    )
    qos, time = ("debug", "00:30:00") if test else ("high", "12:00:00")
    write_job_script(run_folder / "job.sbatch", run_folder, run_folder / "input.txt",
                     threads=threads, qos=qos, time=time, memory="32G",
                     job_name=f"{site.name}_{day:%Y%m%d}", before_name=write_main,
                     name_dir=PATCHED_NAME_DIR)  # fmt: skip
    timing = run_timing(plan, settings)  # unaffected by the pressures
    row.update(
        status="prepared",
        release_groups=members["release_group"].nunique(),
        sources=len(plan),
        met_version=settings.met_version.name,
        met_start=f"{timing['window_start']:%Y-%m-%dT%H:%M}",
        met_end=f"{plan['release_start'].max():%Y-%m-%dT%H:%M}",
    )
    return row


def write_main_input(run_folder):
    """Write a day's release groups, release plan and NAME input from its surface-pressure run.

    Run by the day's main job after the surface-pressure job (``submit.sh``); can be run by hand.

    Writes ``surface_pressure.csv`` (retrieval and model surface pressure, hPa, and their ratio),
    ``release_groups.csv`` (each sounding's release group, D20), ``release_plan.csv`` (one set of
    releases per group, levels scaled to the model surface, D19) and ``input.txt``.

    Args:
        run_folder: the day folder (``WORK/<month folder>/YYYYMMDD``).

    Returns:
        The release plan (DataFrame, as ``release_plan.csv``).

    Raises:
        FileNotFoundError if a surface-pressure output file is missing. ValueError if the model
        pressures are implausible (``surface_pressure.read_surface_pressure``,
        ``surface_pressure_table``), or a modelled
        release layer reaches above the domain top or the modelled levels carry too little of a
        sounding's column weight (``releases.check_modelled_levels``, D21).
    """
    run_folder = Path(run_folder)
    info = json.loads((run_folder / "day.json").read_text())
    site = SITES[info["site"]]
    day = pd.Timestamp(info["day"])
    with xr.open_dataset(run_folder / f"obs_{site.name}_{day:%Y%m%d}.nc") as f:
        obs = f.load()

    folder = run_folder / "surface_pressure"
    soundings = pd.read_csv(folder / "soundings.csv", parse_dates=["time"],
                            dtype={"sounding_id": str})  # fmt: skip
    model = read_surface_pressure(folder / "output", soundings)
    surface = surface_pressure_table(obs, model)
    surface.to_csv(run_folder / "surface_pressure.csv", index=False)

    levels = obs["pressure_levels"].values
    scaled = np.array([scale_to_model_surface(levels[p], surface["model_hpa"].iloc[p])
                       for p in range(len(surface))])  # fmt: skip
    n_levels = MODELLED_LEVELS[info["product"]]
    check_modelled_levels(obs, scaled, surface["domain_top_hpa"].to_numpy(), n_levels)
    rates = PARTICLE_RATES[info["product"]]
    plan, members = grouped_release_plan(obs, site, n_levels, pressure_levels_hpa=scaled,
                                         particle_rates=rates)  # fmt: skip
    members.to_csv(run_folder / "release_groups.csv", index=False)
    plan.to_csv(run_folder / "release_plan.csv", index=False)
    settings = run_settings(site, run_folder, day, info["threads"], info["test"])
    write_name_input(run_folder / "input.txt", plan, settings)
    print(f"{run_folder}: {len(surface)} soundings in {members['release_group'].nunique()} release "
          f"groups, {len(plan)} releases; model/retrieval surface pressure "
          f"{surface['ratio'].min():.4f}-{surface['ratio'].max():.4f}")  # fmt: skip
    return plan


def prepare_month(product, site_name, month, work, threads=THREADS, test=False, days=None):
    """Prepare every day of ``month`` ("YYYY-MM"), or only ``days`` (dates in the month).

    Creates ``work/<site>_<product>_<YYYYMM>/`` with the day folders, ``summary.csv``,
    ``restore_met.sbatch`` and ``submit.sh`` (see the module docstring). Submits nothing.

    Args:
        product: L2 product name, e.g. "GOSAT-1 v10.0".
        site_name: key of ``domains.SITES``, e.g. "GOSAT-BRAZIL".
        month: "YYYY-MM".
        work: parent folder for the month folder (scratch).
        threads, test: as ``prepare_day``.
        days: optional list of dates in the month.

    Returns:
        The summary table (DataFrame, as summary.csv).

    Raises:
        SystemExit if the month folder already exists or no day has soundings; ValueError if a day
        is outside the month or the prepared days need more than one met version.
    """
    site = SITES[site_name]
    month = pd.Period(month, freq="M")
    month_folder = Path(work) / f"{site.name}_{product.replace(' ', '_')}_{month:%Y%m}"
    if month_folder.exists():
        raise SystemExit(f"{month_folder} already exists; choose a new folder")
    month_folder.mkdir(parents=True)

    all_days = pd.date_range(month.start_time, month.end_time.normalize(), freq="D")
    if days is not None:
        chosen = pd.DatetimeIndex(pd.to_datetime(list(days))).normalize()
        if not chosen.isin(all_days).all():
            raise ValueError(f"Days {list(days)} are not all in {month}")
        all_days = chosen
    rows = [prepare_day(product, site, day, month_folder, threads, test) for day in all_days]
    summary = pd.DataFrame(rows)
    summary.to_csv(month_folder / "summary.csv", index=False)
    prepared = summary[summary["status"] == "prepared"]
    if prepared.empty:
        print(summary.to_string(index=False))
        raise SystemExit("No day had soundings to run")

    versions = prepared["met_version"].unique()
    if len(versions) != 1:
        raise ValueError(f"The month spans met versions {list(versions)}: not supported yet")
    met_folder = MET_ROOT / versions[0]
    write_met_restore_job(
        month_folder / "restore_met.sbatch",
        month_folder,
        versions[0],
        prepared["met_start"].min(),
        prepared["met_end"].max(),
        met_folder,
    )
    lines = [
        "#!/bin/bash",
        "# Restore the met; then for each day the surface-pressure run, then the main run.",
        "# afterok: a job starts only if the one before it succeeded.",
        "set -e",
        f"cd {month_folder}",
        "restore=$(sbatch --parsable restore_met.sbatch)",
        'echo "met restore: job $restore"',
    ]
    for day in prepared["day"]:
        folder = pd.Timestamp(day).strftime("%Y%m%d")
        lines += [
            f"ps=$(sbatch --parsable --dependency=afterok:$restore "
            f"{folder}/surface_pressure/job.sbatch)",
            f"main=$(sbatch --parsable --dependency=afterok:$ps {folder}/job.sbatch)",
            f'echo "{folder}: surface pressure job $ps, main job $main"',
        ]
    (month_folder / "submit.sh").write_text("\n".join(lines) + "\n")
    print(summary.to_string(index=False))
    print(f"Prepared {len(prepared)} days in {month_folder}; submit with: bash submit.sh")
    return summary


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if args and args[0] == "main-input":
        write_main_input(args[1])
    else:
        days = [a.removeprefix("--days=") for a in sys.argv[1:] if a.startswith("--days=")]
        prepare_month(*args[:4], test="--test" in sys.argv[1:],
                      days=days[0].split(",") if days else None)  # fmt: skip
