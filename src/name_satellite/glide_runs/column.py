"""The Lamont column release (GLIDE spec §5): model surface pressure each hour, then the plan.

Where this sits
---------------
    **1. met-only NAME runs: model surface pressure at Lamont at every release start**
    ->  **2. column release plan**  ->  prepare.py (NAME input, job)  ->  LOTUS2

Spec §5: every Lamont release is a column, with particles starting at pressures drawn evenly
between the MODEL's surface pressure at the release's start time and 100 hPa, 100,800 particles
(28/s for an hour). A NAME source in pressure (``Cuboid``, ``Z-Coord Pa``) spreads its particles
evenly in pressure (test of 16 Nov 2026, job 2296900; dev/notes/name-behaviour.md), so each release
is ONE source, not the 28 stacked layers the spec allows as a fallback.

1. Surface pressure. ``write_surface_pressure_runs`` writes short met-only NAME runs (one dummy
   particle, ``surface_pressure.surface_pressure_input_text``) that report "Pressure (Pa)" at
   0 m agl, and at the domain top, at the receptor at each release start time. NAME interpolates
   the UM fields to the point and time. One run handles at most ``HOURS_PER_RUN`` times (NAME
   allows 200 time grids), so a 336-hour block needs two. Files, per run folder
   ``<folder>/surface_pressure_<k>/``: ``soundings.csv`` (the times sampled), ``input.txt``,
   ``job.sbatch``, ``restore_met.sbatch`` and ``output/ReleaseMet_S<nnnn>_C1.txt``.
   ``read_surface_pressure_runs`` reads them all into one table (``surface_pressure.csv``).
2. Plan. ``column_release_plan`` makes one source per release from that table.

Command line (from the repository folder)::

    uv run python -m name_satellite.glide_runs.column write FOLDER SGP FIRST LAST MET_FOLDER
    # submit each surface_pressure_<k>/restore_met.sbatch, then its job.sbatch (afterok)
    uv run python -m name_satellite.glide_runs.column read FOLDER   # -> FOLDER/surface_pressure.csv

FIRST and LAST are the first and last release start times (UTC, e.g. 2019-07-01T00).
"""

import sys
from pathlib import Path

import pandas as pd

from name_satellite import lotus
from name_satellite.glide_runs import spec
from name_satellite.release_met import MAX_PER_GROUP
from name_satellite.surface_pressure import (
    TRAVEL_HOURS,
    read_surface_pressure,
    surface_pressure_input_text,
)

HOURS_PER_RUN = 168  # below NAME's 200 time grids per run (release_met.MAX_SOUNDINGS = 199)


def sample_times_table(receptor, release_times):
    """Where and when to sample the surface pressure: one row per release start.

    Returns a ``release_met``-style table: ``sounding`` (S0001, ...), ``release_time``, ``lon``,
    ``lat``, ``time`` (= the release start, §5) and ``group`` (output groups of at most
    ``release_met.MAX_PER_GROUP`` points).
    """
    times = pd.DatetimeIndex(release_times)
    table = pd.DataFrame(
        {
            "sounding": [f"S{i:04d}" for i in range(1, len(times) + 1)],
            "release_time": times,
            "lon": receptor.lon,
            "lat": receptor.lat,
            "time": times,
        }
    )
    table["group"] = table.index // MAX_PER_GROUP + 1
    return table


def write_surface_pressure_runs(folder, receptor_code, release_times, met_folder):
    """Write the met-only runs for ``release_times`` at a column receptor. Returns their folders.

    Args:
        folder: parent folder (each run goes in ``surface_pressure_<k>``; they must not exist).
        receptor_code: a column receptor of ``spec.RECEPTORS`` (Lamont, "SGP").
        release_times: release start times (UTC).
        met_folder: the met folder (shared with the main runs of the block).
    """
    receptor = spec.RECEPTORS[receptor_code]
    if not receptor.is_column:
        raise ValueError(f"{receptor_code} is not a column receptor")
    times = pd.DatetimeIndex(release_times)
    folders = []
    for k, start in enumerate(range(0, len(times), HOURS_PER_RUN), start=1):
        run = Path(folder) / f"surface_pressure_{k}"
        run.mkdir(parents=True)
        (run / "output").mkdir()
        table = sample_times_table(receptor, times[start : start + HOURS_PER_RUN])
        table.to_csv(run / "soundings.csv", index=False)
        settings = spec.reference_settings(
            run_name=f"{receptor_code}_surface_pressure_{k}",
            output_folder=run / "output",
            domain=spec.DOMAINS[receptor.domain],
            met_folder=met_folder,
            max_particles=100,
            threads=1,
        )
        (run / "input.txt").write_text(surface_pressure_input_text(table, settings))
        lotus.write_job_script(
            run / "job.sbatch", run, run / "input.txt", threads=1, qos="short", time="01:00:00",
            memory="8G", job_name=f"glide_{receptor_code}_ps{k}",
        )  # fmt: skip
        # Met from the first sample (less the dummy particle's travel) to the last, ±1 h.
        lotus.write_met_restore_job(
            run / "restore_met.sbatch",
            run,
            settings.met_version.name,
            (table["time"].min() - pd.Timedelta(hours=TRAVEL_HOURS + 1)).strftime("%Y-%m-%dT%H:%M"),
            (table["time"].max() + pd.Timedelta(hours=1)).strftime("%Y-%m-%dT%H:%M"),
            met_folder,
            time="02:00:00",
        )
        folders.append(run)
    return folders


def read_surface_pressure_runs(folders):
    """Model surface pressure at every sampled release start, from the met-only runs.

    Returns:
        DataFrame, one row per release, sorted by time: ``release_time``,
        ``surface_pressure_pa`` and ``domain_top_pressure_pa`` (model pressure at the domain's
        top, 20 km agl).

    Raises:
        ValueError if a release time appears twice, or the column top (100 hPa) is not inside the
        domain (its pressure must be above the domain top's).
    """
    parts = []
    for run in folders:
        table = pd.read_csv(Path(run) / "soundings.csv", parse_dates=["release_time", "time"])
        parts.append(read_surface_pressure(Path(run) / "output", table))
    table = pd.concat(parts, ignore_index=True).sort_values("release_time")
    if table["release_time"].duplicated().any():
        raise ValueError("A release time was sampled twice")
    top_pa = spec.COLUMN_TOP_HPA * 100.0
    if not (table["domain_top_pressure_pa"] < top_pa).all():
        raise ValueError("The column top (100 hPa) is above the domain top at some time")
    columns = ["release_time", "surface_pressure_pa", "domain_top_pressure_pa"]
    return table[columns].reset_index(drop=True)


def column_release_plan(receptor, surface_pressure):
    """Release plan for the column receptor: one source per release, even in pressure (§5).

    Args:
        receptor: the column ``spec.Receptor`` (Lamont).
        surface_pressure: ``read_surface_pressure_runs`` table (one row per release).

    Returns:
        DataFrame in the columns of ``plan.point_release_plan``, with ``z_coord`` "Pa", ``z`` the
        middle of the column and ``dz`` its full depth (Pa), plus ``pressure_bottom_pa`` and
        ``pressure_top_pa`` for the delivered file (§8, in hPa there).
    """
    if not receptor.is_column:
        raise ValueError(f"{receptor.code} is not a column receptor")
    times = pd.DatetimeIndex(surface_pressure["release_time"])
    if (times != times.floor("1h")).any() or not times.is_monotonic_increasing:
        raise ValueError("Release times must be on the hour and increasing")
    bottom = surface_pressure["surface_pressure_pa"].to_numpy(dtype=float)
    top = spec.COLUMN_TOP_HPA * 100.0
    if not (bottom > top).all():
        raise ValueError("Surface pressure below the column top")
    hour = pd.Timedelta(hours=spec.RELEASE_HOURS)
    n = len(times)
    return pd.DataFrame(
        {
            "source_name": [f"SourceID_{i}" for i in range(1, n + 1)],
            "release_time": times,
            "release_start": times + hour,
            "release_end": times,
            "receptor": receptor.code,
            "lon": receptor.lon,
            "lat": receptor.lat,
            "z_coord": "Pa",
            "z": (bottom + top) / 2.0,
            "dz": bottom - top,
            "dlon": 0.0,
            "dlat": 0.0,
            "particles_per_second": spec.COLUMN_PARTICLES_PER_SECOND,
            "pressure_bottom_pa": bottom,
            "pressure_top_pa": top,
        }
    )


# For reference in the run record: how the column is built.
COLUMN_METHOD = (
    "one NAME Cuboid source per release, Z-Coord Pa, from the model surface pressure at the "
    "release start (NAME point met, interpolated to Lamont) to 100 hPa; NAME spreads particles "
    "evenly in pressure (measured 16 Nov 2026, job 2296900)"
)


def main(argv):
    """``write FOLDER RECEPTOR FIRST LAST MET_FOLDER`` or ``read FOLDER`` (see the docstring)."""
    action, folder = argv[0], Path(argv[1])
    if action == "write":
        receptor, first, last, met_folder = argv[2:6]
        times = pd.date_range(first, last, freq="1h")
        for run in write_surface_pressure_runs(folder, receptor, times, met_folder):
            print(f"Wrote {run}")
    elif action == "read":
        runs = sorted(folder.glob("surface_pressure_*/"))
        table = read_surface_pressure_runs(runs)
        table.to_csv(folder / "surface_pressure.csv", index=False)
        print(f"Wrote {folder / 'surface_pressure.csv'}: {len(table)} release times")
    else:
        raise SystemExit("Usage: column write FOLDER RECEPTOR FIRST LAST MET_FOLDER | read FOLDER")


if __name__ == "__main__":
    main(sys.argv[1:])
