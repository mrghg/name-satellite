"""Run the LEGACY NAME_processing `process()` on a copied workspace.

Purpose: check whether today's legacy processing code reproduces the processed footprint files
that were made in 2018 (concentration run) and 2019 (mixing-ratio run). If it does, those files
become exact month-scale tests for our new processing code.

This script uses the legacy code, so it runs in the legacy environment, not in uv:

    module load jaspy/3.12/v20250704
    PYTHONPATH=~/code/NAME_processing python dev/scripts/legacy_validation/run_legacy_process.py WORKSPACE FIELDS_FOLDER

WORKSPACE must be a COPY of the inputs in scratch (never `data/`), laid out as the legacy code
expects:
    WORKSPACE/SOUTHAMERICA_GOSAT-BRAZIL_column/{FIELDS_FOLDER,Particle_files,Met,Observations}
    WORKSPACE/SOUTHAMERICA_GOSAT-BRAZIL_column/time_step.txt   (release length in seconds: 60)
FIELDS_FOLDER is "MixR_files" (mixing-ratio run) or "Fields_files" (concentration run).

The processed file is written to WORKSPACE/processed/.
The settings match the legacy test suite (NAME_processing/tests/test_process_sat.py:
upper_level = max_level = 17).
"""

import sys
from pathlib import Path

import name_processing.process as legacy

workspace = Path(sys.argv[1]).resolve()
fields_folder = sys.argv[2]

if "/code/name-satellite/data" in str(workspace):
    raise SystemExit("Refusing to run on the reference data folder; use a scratch copy.")

output_dir = workspace / "processed"
output_dir.mkdir(exist_ok=True)

result = legacy.process(
    domain="SOUTHAMERICA",
    site="GOSAT-BRAZIL",
    height="column",
    year=2013,
    month=1,
    data_type="satellite",
    upper_level=17,
    max_level=17,
    base_dir=workspace,
    process_dir=output_dir,
    fields_folder=fields_folder,
    particles_folder="Particle_files",
    met_folder="Met",
    processed_folder=str(output_dir),
    obs_folder="Observations",
)
print("process() returned:", type(result))
print("Output files:", sorted(str(p) for p in output_dir.rglob("*.nc")))
