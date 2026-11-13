"""Legacy processing of every January 2013 day of the 2019 MR run, compared with production.

Runs ``legacy_pipeline.process_legacy_inputs`` on each day's raw output of the 2019 mixing-ratio
run (the run the production footprints came from), writes the result with
``write_legacy_fp_file``, and checks every variable against the production file
(``refdata`` "production_201301") for exact equality. The release length is 60 s (the
archive has no time_step.txt).

    uv run python dev/scripts/legacy_validation/check_month_against_production.py RUN_FOLDER OUTPUT_FOLDER

RUN_FOLDER is the ``SOUTHAMERICA_GOSAT-BRAZIL_column`` folder (MixR_files, Particle_files, Met,
Observations). OUTPUT_FOLDER must not exist yet.
"""

import sys
from pathlib import Path

import numpy as np
import xarray as xr

from name_satellite.name_output import read_acrg_particle_file, read_fields
from name_satellite.validation.legacy_files import read_legacy_met_folder
from name_satellite.validation.legacy_pipeline import process_legacy_inputs, write_legacy_fp_file
from name_satellite.validation.refdata import ref_path

TAG = "SOUTHAMERICA_GOSAT-BRAZIL_column"
RELEASE_SECONDS = 60.0


def check_day(run, out, day, production):
    """Process one day; return a one-line verdict comparing every variable with production."""
    fields = read_fields(run / "MixR_files" / f"Fields_{TAG}_{day}.txt.gz")
    particles = read_acrg_particle_file(
        run / "Particle_files" / f"particle_location_{TAG}_{day}.txt.gz"
    )
    met = read_legacy_met_folder(run / "Met" / f"Met_{TAG}_{day}")
    obs = xr.open_dataset(run / "Observations" / f"gosat-fts_gosat_{day}_ch4-column.nc").load()
    result = process_legacy_inputs(fields, particles, met, obs, release_seconds=RELEASE_SECONDS)

    path = out / f"fp_{day}.nc"
    write_legacy_fp_file(result, path)
    mine = xr.open_dataset(path).load()
    prod = production.sel(time=mine["time"].values)

    different = []
    for name in mine.data_vars:
        a = mine[name].values
        b = prod[name].transpose(*mine[name].dims).values
        if not np.array_equal(a, b, equal_nan=True):
            relative = np.nanmax(np.abs(a - b)) / np.nanmax(np.abs(b))
            different.append(f"{name} (max diff / max {relative:.1e})")
    verdict = "IDENTICAL" if not different else "DIFFER: " + ", ".join(different)
    return f"{day}: {mine.sizes['time']} soundings, {verdict}"


def main(run, out):
    run, out = Path(run), Path(out)
    if out.exists():
        raise SystemExit(f"{out} already exists; choose a new folder")
    out.mkdir(parents=True)
    production = xr.open_dataset(ref_path("production_201301")).load()
    for d in range(1, 32):
        print(check_day(run, out, f"201301{d:02d}", production), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
