"""Run the LEGACY GOSAT preprocessing (acrg.satellite.gosat) on one L2 file, as a reference.

Our observation reader (name_satellite.observations) is tested against the output of this script.
The settings are those recorded in the 2013 legacy observation file (docs/reference-data.md):
quality filter, bad-pressure filter, binning to 0.234 deg × 0.352 deg within fixed lat/lon bounds.
The NAME-surface-pressure filter is NOT applied (it needs legacy surface-pressure files).

The legacy code runs in the legacy environment, with a few packages jaspy lacks installed into a
scratch folder (uv pip install --target, never into the system):

    module load jaspy/3.12/v20250704
    PYTHONPATH=/work/scratch-pw5/mattrigby/name-satellite/legacy_pydeps_20261109:~/code/acrg \
        python dev/scripts/run_legacy_gosat.py L2_FILE OUTPUT_FOLDER

Writes, under OUTPUT_FOLDER (which must not exist yet), the dataset that gosat_process_file
returns after filtering and binning (before its own file writer, which fails with today's xarray
because it expects the exposure IDs as str, not bytes):
    binned.nc     filtered and binned
    unbinned.nc   filtered, lat/lon bounds only
Python's random module is seeded, so the legacy random time offsets for duplicate times (X12) are
repeatable. The legacy code cannot process GOSAT-1 v10 or GOSAT-2 v2 files: their times are whole
minutes, so many soundings share a time and the binning fails (docs/deviations.md X22).
"""

import random
import sys
from pathlib import Path

import acrg.satellite.gosat as legacy

LAT_BOUNDS = [-35.753, 7.251]  # from the 2013 legacy obs file's "bins" attribute
LON_BOUNDS = [-75.984, -32.793]
COORD_BIN = [0.234, 0.352]  # [lat, lon] degrees

l2_file = Path(sys.argv[1]).resolve()
out = Path(sys.argv[2]).resolve()
if "/code/name-satellite/data" in str(out):
    raise SystemExit("Refusing to write into the reference data folder")
if out.exists():
    raise SystemExit(f"{out} already exists; choose a new folder")

out.mkdir(parents=True)
for name, coord_bin in (("binned", COORD_BIN), ("unbinned", None)):
    random.seed(0)
    result = legacy.gosat_process_file(
        str(l2_file),
        site="GOSAT-BRAZIL",
        species="ch4",
        lat_bounds=LAT_BOUNDS,
        lon_bounds=LON_BOUNDS,
        coord_bin=coord_bin,
        quality_filt=True,
        bad_pressure_filt=True,
        name_sp_filt=False,
        mode=None,
        use_name_pressure=False,
        write_nc=False,
        write_name=False,
        verbose=False,
    )
    # Exposure IDs: bytes (unbinned) or comma-joined text (binned); store as text.
    ids = [x.decode() if isinstance(x, bytes) else str(x) for x in result["exposure_id"].values]
    result["exposure_id"] = ("time", ids)
    result.to_netcdf(out / f"{name}.nc")
    print(f"{name}: {result.sizes['time']} soundings")
