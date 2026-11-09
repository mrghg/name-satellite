"""Prepare a tiny NAME 8.7.1 test run on LOTUS2, to check that our input files work.

Two releases (sounding 1 of 1 Jan 2013, levels 1-2), 3 h maximum particle age, legacy surface
pressure (X9), Mk6 met. Writes into a NEW folder (refuses to overwrite):

    <work>/input.txt       NAME input file
    <work>/output/         NAME output folder
    <work>/job.sbatch      Slurm job script (debug queue, account acrg)

Met is restored by NAME's restore script into a shared folder <scratch>/met/UMG_Mk6PT.

    uv run python dev/scripts/prepare_debug_run.py /work/scratch-pw5/mattrigby/name-satellite/debug_20261109
    sbatch /work/scratch-pw5/mattrigby/name-satellite/debug_20261109/job.sbatch
"""

import sys
from pathlib import Path

import xarray as xr

from name_satellite.legacy import legacy_surface_pressure_hpa, read_name_input_sources
from name_satellite.met import NAME_DIR, met_version_for_period
from name_satellite.name_input import DispersionDomain, OutputGrid, RunSettings, write_name_input
from name_satellite.refdata import ref_path
from name_satellite.releases import build_release_plan

SCRATCH = Path("/work/scratch-pw5/mattrigby/name-satellite")

JOB_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=name_debug
#SBATCH --account=acrg
#SBATCH --partition=debug
#SBATCH --qos=debug
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output={work}/slurm-%j.out
#SBATCH --error={work}/slurm-%j.err

# Environment as in the NAME team's 8.7.1 demo scripts
module load jaspy
module load oneapi/compilers/24.2.0
module load oneapi/mpi/24.2.0
module load netcdf/intel2024.2.0/4.9.2
module load netcdf/intel2024.2.0/fortran/4.6.1
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_STACKSIZE='32m'
# Shumlib (needed to read UM met files), as in the demo scripts
export LD_LIBRARY_PATH=${{LD_LIBRARY_PATH}}:/gws/ssde/j25a/name/lib/shumlib_2024.11.1/jasmin-lotus-ifx-gcc/lib

cd {work}
echo "Started $(date)"
time mpirun -np 1 {name_dir}/Executables_Linux/nameiii_64bit_par.exe {work}/input.txt
echo "NAME exit code $?; finished $(date)"
"""


def main(work):
    work = Path(work)
    if work.exists():
        raise SystemExit(f"{work} already exists; choose a new folder")
    (work / "output").mkdir(parents=True)
    met_folder = SCRATCH / "met" / "UMG_Mk6PT"
    met_folder.mkdir(parents=True, exist_ok=True)

    obs = xr.open_dataset(ref_path("obs"))
    legacy_sources = read_name_input_sources(ref_path("mr2019_input"))
    plan = build_release_plan(obs, surface_pressure_hpa=legacy_surface_pressure_hpa(legacy_sources))
    plan = plan.iloc[:2]  # sounding 1, levels 1 and 2

    settings = RunSettings(
        run_name="debug_20130101",
        output_folder=str(work / "output"),
        grid=OutputGrid(nx=190, ny=357, dx=0.352, dy=0.234, x_min=-91.329, y_min=-60.981),
        domain=DispersionDomain(x_min=-96.329, x_max=-19.449, y_min=-64.981, y_max=27.557),
        met_version=met_version_for_period("2013-01-01", "2013-01-01"),
        met_folder=str(met_folder),
        max_age_hours=3,
    )
    write_name_input(work / "input.txt", plan, settings)
    (work / "job.sbatch").write_text(JOB_TEMPLATE.format(work=work, name_dir=NAME_DIR))
    print(f"Prepared {work}")


if __name__ == "__main__":
    main(sys.argv[1])
