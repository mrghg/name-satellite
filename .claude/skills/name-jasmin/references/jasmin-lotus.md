# Running NAME on JASMIN / LOTUS2

Working job-script writer: `~/code/name-satellite/src/name_satellite/lotus.py`.

## Environment (as the NAME team's 8.7.1 demo scripts)
```bash
module load jaspy
module load oneapi/compilers/24.2.0
module load oneapi/mpi/24.2.0
module load netcdf/intel2024.2.0/4.9.2
module load netcdf/intel2024.2.0/fortran/4.6.1
export OMP_NUM_THREADS=4          # = --cpus-per-task = input Threads
export OMP_PROC_BIND=close
export OMP_STACKSIZE='32m'
ulimit -s unlimited               # docs recommend a large stack on Linux (A01 §1.8.4); our runs worked without
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/gws/ssde/j25a/name/lib/shumlib_2024.11.1/jasmin-lotus-ifx-gcc/lib
mpirun -np 1 /gws/ssde/j25a/name/code/NAME-8.7.1/Executables_Linux/nameiii_64bit_par.exe input.txt
```
Without the shumlib path NAME stops at start-up: `libshum_wgdos_packing.so: cannot open shared
object file`. Batch jobs don't load your shell's PATH: call tools by full path (e.g.
`~/.local/bin/uv`).

## Slurm
Every job needs `--account`, `--partition` and `--qos`; a wrong pairing fails with "Invalid
account or account/partition combination specified". Submit from the sci servers.

| QoS | Partition | Max CPUs | Max wall | Use |
|---|---|---|---|---|
| `debug` | `debug` | 8 | 1 h | tests (1-h travel runs) |
| `short` | `standard` | 1 | 4 h | single CPU |
| `standard` | `standard` | 1 | 24 h | single CPU |
| `long` | `standard` | 1 | 5 days | single CPU |
| `high` | `standard` | 96 | 48 h | multi-threaded NAME, big memory |

- 30-day back-run of 352 releases on 4 threads: ~82 min, 5.4 GB; ask `--mem=32G`, `--time=12:00:00`.
- Chain jobs with `--dependency=afterok:<id>` (met restore → NAME runs); with a job script that
  fails on a non-empty error file, a bad day stops its dependants.
- Job scripts that worked: `#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=N --hint=compute_bound`,
  output/error to `<work>/slurm-%j.out/.err`, `cd` into the run folder, then mpirun, then check
  `<input>Error.txt`.
- Legacy ACRG job scripts use retired settings (`par-single`, `ivybridge`, `intel/17`,
  `/gws/nopw/j04` paths) and don't work.
- Waiting: check `squeue -u $USER` / `sacct -j <ids> -X -o JobID,JobName,Elapsed,State` now and
  then; if looping, at most once a minute.

## Measuring NAME's peak memory
`sacct` MaxRSS misses short jobs, and `/usr/bin/time -v mpirun ...` measures mpirun (10 MB), not
NAME. What works: run mpirun in the background and read the peak resident memory (`VmHWM` in
`/proc/<pid>/status`) of processes whose command line contains this run's input path, every 2 s,
writing the maximum to a file. Use the input path to tell runs apart when several share a node.

## Storage
| Where | Notes |
|---|---|
| Home `~` | 100 GB quota: code and small files only |
| ACRG GWS `/gws/ssde/j25b/acrg/<user>` | nearly full: keep only processed results |
| Scratch `/work/scratch-pw5/<user>` | NAME runs, restored met; files deleted 28 days after last access |

Avoid more than ~100,000 files in one directory. One day's 30-day back-run on Mk6 needs ~3,500
met files (14 parts × 3-hourly × 31 days, plus one step each end); a month shares most of them.

## Access
Unix groups `gws_name` (NAME workspace, met, patched build) and `gws_acrg`; check with `groups`.
