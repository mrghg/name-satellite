"""Step 4 of the workflow: Slurm job scripts for running NAME on LOTUS2 (JASMIN).

Where this sits (see docs/workflow.md)
--------------------------------------
    3. write NAME input  ->  **4. run NAME on LOTUS2**  ->  5. read NAME output

The environment (modules, OpenMP variables, shumlib library path) is the one used by the NAME
team's 8.7.1 demo scripts; see docs/jasmin.md.

Threads: NAME will not use more threads than OpenMP reports as available
(``OMP_NUM_THREADS``), and LOTUS2 only gives a job the CPUs it asks for
(``--cpus-per-task``). So the job's CPUs, ``OMP_NUM_THREADS`` and the input file's ``Threads``
must all be the same number; here they all come from ``RunSettings.threads``.
"""

from pathlib import Path

from name_satellite.met import NAME_DIR, NAME_GWS

# Maximum CPUs per job for each LOTUS2 QoS (help.jasmin.ac.uk, checked 2026-11-09), and the
# partition each QoS belongs to.
QOS_LIMITS = {
    "debug": {"partition": "debug", "max_cpus": 8},
    "short": {"partition": "standard", "max_cpus": 1},
    "standard": {"partition": "standard", "max_cpus": 1},
    "long": {"partition": "standard", "max_cpus": 1},
    "high": {"partition": "standard", "max_cpus": 96},
}

SHUMLIB = NAME_GWS / "lib" / "shumlib_2024.11.1" / "jasmin-lotus-ifx-gcc" / "lib"

JOB_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --account={account}
#SBATCH --partition={partition}
#SBATCH --qos={qos}
#SBATCH --time={time}
#SBATCH --mem={memory}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={threads}
#SBATCH --hint=compute_bound
#SBATCH --output={work}/slurm-%j.out
#SBATCH --error={work}/slurm-%j.err

# NAME 8.7.1 environment, as in the NAME team's demo scripts (docs/jasmin.md)
module load jaspy
module load oneapi/compilers/24.2.0
module load oneapi/mpi/24.2.0
module load netcdf/intel2024.2.0/4.9.2
module load netcdf/intel2024.2.0/fortran/4.6.1
export OMP_NUM_THREADS={threads}  # must equal --cpus-per-task and the input file's Threads
export OMP_PROC_BIND=close
export OMP_STACKSIZE='32m'
export LD_LIBRARY_PATH=${{LD_LIBRARY_PATH}}:{shumlib}  # needed to read UM met files
{before_name}
cd {work}
echo "Started $(date) on $(hostname) with {threads} thread(s)"
time mpirun -np 1 {name_dir}/Executables_Linux/nameiii_64bit_par.exe {input_file}
status=$?
echo "NAME exit code $status; finished $(date)"
# NAME can exit 0 although it reported errors (e.g. lost particles): judge the run by its error
# file too (docs/jasmin.md). A non-zero exit makes dependent jobs (afterok) not start.
error_file={error_file}
if [ -s "$error_file" ]; then
    echo "Error file $error_file is not empty:"; cat "$error_file"
    [ $status -eq 0 ] && status=1
fi
exit $status
"""

# Restoring (copying and unzipping) the met a run needs from the NAME met archive, as its own job:
# NAME jobs must not restore into a shared met folder at the same time (dev/notes, benchmark).
MET_RESTORE_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --account={account}
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time={time}
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output={work}/restore-%j.out
#SBATCH --error={work}/restore-%j.err

cd {repo}
time ~/.local/bin/uv run --no-sync python -m name_satellite.met_restore \\
    {version} {start} {end} {folder}
"""

REPO = Path(__file__).resolve().parents[2]


def write_job_script(
    path,
    work_dir,
    input_file,
    threads,
    qos="debug",
    time="00:30:00",
    memory="16G",
    job_name="name_satellite",
    account="acrg",
    before_name="",
):
    """Write a Slurm job script that runs NAME 8.7.1 on one input file.

    Args:
        path: job script to write.
        work_dir: folder the job runs in (Slurm logs go here).
        input_file: NAME input file (written with the same number of threads).
        threads: OpenMP threads = CPUs requested.
        qos: LOTUS2 QoS; the partition follows from it. Must allow ``threads`` CPUs.
        time, memory: Slurm limits, e.g. "02:00:00", "16G".
        account: Slurm account (decision D15: acrg).
        before_name: shell lines to run before NAME (e.g. writing the input file,
            ``prepare.write_main_input``); the job stops with an error if they fail.

    Returns:
        The path written.
    """
    if qos not in QOS_LIMITS:
        raise ValueError(f"Unknown QoS {qos!r}; choose from {sorted(QOS_LIMITS)}")
    if not 1 <= threads <= QOS_LIMITS[qos]["max_cpus"]:
        raise ValueError(
            f"QoS {qos!r} allows at most {QOS_LIMITS[qos]['max_cpus']} CPU(s); asked for {threads}"
        )
    text = JOB_TEMPLATE.format(
        job_name=job_name,
        account=account,
        partition=QOS_LIMITS[qos]["partition"],
        qos=qos,
        time=time,
        memory=memory,
        threads=threads,
        work=work_dir,
        shumlib=SHUMLIB,
        name_dir=NAME_DIR,
        input_file=input_file,
        before_name=f"\n{before_name.rstrip()}\n" if before_name else "",
        # NAME names its error file after the input file: input.txt -> inputError.txt
        error_file=Path(work_dir) / (Path(input_file).stem + "Error.txt"),
    )
    path = Path(path)
    path.write_text(text)
    return path


def write_met_restore_job(
    path, work_dir, met_version_name, start, end, met_folder, time="08:00:00", account="acrg"
):
    """Write a Slurm job that restores the met for a period into a met folder.

    Args:
        path: job script to write.
        work_dir: folder for the job's logs.
        met_version_name: e.g. "UMG_Mk6PT" (``met.MET_VERSIONS``).
        start, end: period of met needed (ISO strings, e.g. "2012-12-02T14:58").
        met_folder: where the met files go (shared by the runs that use it).
        time: Slurm time limit.
        account: Slurm account (D15).

    Returns:
        The path written.
    """
    text = MET_RESTORE_TEMPLATE.format(
        job_name=f"restore_{met_version_name}",
        account=account,
        time=time,
        work=work_dir,
        repo=REPO,
        version=met_version_name,
        start=start,
        end=end,
        folder=met_folder,
    )
    path = Path(path)
    path.write_text(text)
    return path
