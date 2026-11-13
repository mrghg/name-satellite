"""Step 1b of the workflow: the model's surface pressure at each sounding, from a met-only NAME run.

Where this sits (see docs/workflow.md)
--------------------------------------
    1. read soundings  ->  **1b. met-only NAME run: model surface pressure**  ->  2. plan
    releases (scaled to the model surface, D19)  ->  3. write NAME input  ->  4. run NAME

Why
---
The release layers follow the met model's terrain (decision D19, ``releases.py``): each
sounding's retrieval levels are scaled by model / retrieval surface pressure. So the model's
surface pressure at each sounding must be known BEFORE the main NAME input can be written. The
main run's own release-point met (``release_met.py``) comes too late, so this short run gets it.

The run
-------
NAME needs at least one source, so the run has one dummy source of a single particle, 10 m above
the ground at the first sounding, released over the day's soundings and stopped after 1 hour of
travel (the legacy separate met run did the same, ``Met_satellite_template.txt``). It requests
two quantities at each sounding's position and time (``QUANTITIES``): "Pressure (Pa)" at 0 m agl,
i.e. the model's surface pressure (D19), and "Pressure (Pa)" at the dispersion domain's top (20 km
agl for new runs), to check the modelled release layers lie inside the domain (D21). The requests
are written and read with ``release_met``. Same met version, met folder, domain and grid as the
main run; stock NAME 8.7.1; 1 thread; takes minutes.

Files (in ``<day folder>/surface_pressure/``)
---------------------------------------------
* ``soundings.csv`` (written by ``prepare.prepare_day``): one row per sounding, where and when
  to sample (``release_met.release_met_soundings``): ``sounding`` (S0001, ...), ``point_number``,
  ``sounding_id``, ``lon``, ``lat`` (degrees), ``time`` (UTC, truncated to the second), ``z_pa``,
  ``group``.
* ``input.txt`` (``surface_pressure_input_text``) and ``job.sbatch``.
* ``output/ReleaseMet_<sounding>_C1.txt`` (written by NAME): one NAME point-output file per
  sounding. In the header, the ``<sounding>_PS`` / ``<sounding>_PT`` columns are the requests and
  the "Z = ..." row their heights; the data rows (one per sounding location of its output group)
  hold the values in Pa; only the row of the sounding itself is used.

``prepare.write_main_input`` reads them (``read_surface_pressure``) and writes
``<day folder>/surface_pressure.csv`` (``surface_pressure_table``, hPa).
"""

import dataclasses

import numpy as np
import pandas as pd

from name_satellite.name_input import (
    NAME_TIME_FORMAT,
    dispersion_options_block,
    domain_block,
    files_and_coordinates_block,
    main_options_block,
    met_declarations_block,
    openmp_block,
    restart_and_cases_block,
    species_block,
    timer_block,
)
from name_satellite.release_met import read_release_met, release_met_blocks

# Model surface pressure (D19), and the model pressure at the dispersion domain's top, so the
# main job can check that every modelled release layer lies inside the domain (D21).
QUANTITIES = [("PS", "Pressure (Pa)", "ground"), ("PT", "Pressure (Pa)", "top")]
TRAVEL_HOURS = 1  # the dummy particle only has to exist; it is not used

# Model / retrieval surface pressure outside this range means something is wrong (position, time
# or units), not a topography difference: 1 Jan 2013 ranged 0.992-1.036.
PLAUSIBLE_RATIO = (0.9, 1.1)


def surface_pressure_settings(settings, output_folder):
    """Settings for the met-only run, from the main run's (same met, domain and grid).

    Changes only the run name (``<main run name>_surface_pressure``), the output folder, the
    particle limit (100), the threads (1) and switches off the main run's release-met requests
    (this run makes its own).
    """
    return dataclasses.replace(
        settings,
        run_name=f"{settings.run_name}_surface_pressure",
        output_folder=str(output_folder),
        max_particles=100,
        threads=1,
        release_met=False,
    )


def dummy_source_block(soundings):
    """The ``Sources:`` block: one dummy particle (NAME needs a source to run).

    Released from the latest sounding time + 60 s back to the earliest (a backwards run starts at
    the later time), 10 m above the ground at the first sounding of ``soundings`` (a
    ``release_met_soundings`` table). Its output is not used.
    """
    first = soundings.iloc[0]
    start = soundings["time"].max() + pd.Timedelta(seconds=60)  # backwards: start is later
    stop = soundings["time"].min()
    return (
        "Sources:\n"
        "Name, Start Time, Stop Time, # Particles, Source Strength, Shape, H-Coord, Z-Coord, "
        "Set of Locations, Location, dH-Metres?, dZ-Metres?, X, Y, Z, dX, dY, dZ, Angle, "
        "Top Hat, Plume Rise?, Temperature, Volume Flow Rate, Max Age, No Reflect?\n"
        f"DummySource, {start.strftime(NAME_TIME_FORMAT)}, {stop.strftime(NAME_TIME_FORMAT)}, 1, "
        f"INERT-TRACER 1.0 g/s, Cuboid, Lat-Long, m agl, , , No, No, {float(first.lon)!r}, "
        f"{float(first.lat)!r}, 10.0, 0.0, 0.0, 0.0, 0.0, Yes, No, 0.0, 0.0, infinity, Yes\n"
    )


def surface_pressure_input_text(soundings, settings):
    """Text of the met-only NAME input file.

    Args:
        soundings: one row per sounding (``release_met.release_met_soundings`` of a provisional
            release plan): ``sounding``, ``lon``, ``lat``, ``time``, ``group``.
        settings: from ``surface_pressure_settings``.
    """
    blocks = [
        "! Met-only NAME run: model surface pressure at each sounding (name_satellite."
        "surface_pressure)\n",
        main_options_block(settings, 1),
        restart_and_cases_block(),
        openmp_block(settings),
        timer_block(settings),
        files_and_coordinates_block(settings),
        domain_block(settings, {"max_travel_hours": TRAVEL_HOURS}),
        species_block(),
        dispersion_options_block(settings),
        dummy_source_block(soundings),
        met_declarations_block(settings),
        release_met_blocks(
            soundings, settings.grid.dx, settings.grid.dy, QUANTITIES, settings.domain.z_max_m_agl
        ),
    ]
    return "\n".join(blocks)


def read_surface_pressure(output_folder, soundings):
    """Model surface pressure and domain-top pressure (Pa) at each sounding, from the met-only run.

    Args:
        output_folder: the run's output folder (``ReleaseMet_<sounding>_C1.txt`` files).
        soundings: the run's ``soundings.csv`` table (``sounding``, ``sounding_id``, ...).

    Returns:
        ``soundings`` with ``surface_pressure_pa`` and ``domain_top_pressure_pa`` added (Pa).

    Raises:
        FileNotFoundError if an output file is missing. ValueError if a request is missing from
        a file, a value is not a positive number, or the domain top's pressure is not below the
        surface pressure.
    """
    table = read_release_met(output_folder, soundings, QUANTITIES)
    surface = table["surface_pressure_pa"].to_numpy()
    top = table["domain_top_pressure_pa"].to_numpy()
    if not np.all(np.isfinite(surface) & (surface > 0) & np.isfinite(top) & (top > 0)):
        raise ValueError(f"Bad model pressure in {output_folder}")
    if not np.all(top < surface):
        raise ValueError(f"Domain-top pressure not below the surface pressure in {output_folder}")
    return table


def surface_pressure_table(obs, model):
    """Retrieval and model surface pressure per sounding, matched by sounding ID.

    Args:
        obs: observation Dataset (``observations.to_obs``: dim ``sounding_id``).
        model: ``read_surface_pressure`` output, with ``sounding_id``.

    Returns:
        DataFrame in the order of ``obs``: ``sounding_id``, ``retrieval_hpa`` (level 0),
        ``model_hpa``, ``ratio`` (model / retrieval: the factor applied to every level, D19) and
        ``domain_top_hpa`` (model pressure at the dispersion domain's top, D21).

    Raises:
        ValueError if the soundings differ or a ratio is outside PLAUSIBLE_RATIO.
    """
    ids = [str(i) for i in obs["sounding_id"].values]
    model = model.assign(sounding_id=model["sounding_id"].astype(str)).set_index("sounding_id")
    if sorted(model.index) != sorted(ids):
        raise ValueError("Surface-pressure run and observation file have different soundings")
    table = pd.DataFrame(
        {
            "sounding_id": ids,
            "retrieval_hpa": obs["pressure_levels"].values[:, 0].astype(float),
            "model_hpa": model.loc[ids, "surface_pressure_pa"].to_numpy() / 100.0,
        }
    )
    table["ratio"] = table["model_hpa"] / table["retrieval_hpa"]
    table["domain_top_hpa"] = model.loc[ids, "domain_top_pressure_pa"].to_numpy() / 100.0
    low, high = PLAUSIBLE_RATIO
    bad = table[(table["ratio"] < low) | (table["ratio"] > high)]
    if not bad.empty:
        raise ValueError(f"Implausible model/retrieval surface pressure:\n{bad.to_string()}")
    return table
