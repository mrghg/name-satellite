# Running on JASMIN

## Where to run what
- **Sci / login node:** editing, `uv sync`, unit tests, single-day processing, inspecting output.
  Keep runs short and modest in memory; these nodes are shared.
- **LOTUS2 (SLURM):** all NAME runs and anything heavy (month-scale processing, many days).
  Submit with `sbatch`; never run NAME interactively.

## LOTUS2 notes (checked 2026-11-09 against help.jasmin.ac.uk)
Jobs are submitted with `sbatch` directly from the sci servers. Every job needs `--account`,
`--partition` and `--qos` (and should give `--mem` and `--time`); a wrong pairing fails with
"Invalid account or account/partition combination specified".

- **Account: use `--account=acrg`** (decided 2026-11-09). Matt is also in `name`, `footprints`,
  `inhale`, `dare_uk`, `ai4netzero` (ignore `shobu`). The account only decides which project the
  compute is charged to; file access comes from Unix groups (`gws_name`, `gws_acrg`,
  `gws_footprints`), which every job inherits.
- **Partition / QoS:**

  | QoS | Partition | Max CPUs/job | Max wall time | Max memory | Max jobs/user | Use |
  |---|---|---|---|---|---|---|
  | `debug` | `debug` | 8 | 1 h | – | 32 | testing; dedicated nodes |
  | `short` | `standard` | 1 | 4 h | 200 GB | 2000 | short single-CPU jobs |
  | `standard` | `standard` | 1 | 24 h | 128 GB | 4000 | default |
  | `long` | `standard` | 1 | 5 days | 128 GB | 1350 | long, lower priority |
  | `high` | `standard` | 96 | 48 h | 1000 GB | – | multi-CPU / big memory |

  Note: only `debug` and `high` allow more than one CPU, so multi-threaded NAME runs need `high`
  (or `debug` for tests, up to 8 CPUs).
- Job arrays: at most 10,000 elements. Do **not** poll Slurm with `watch` or similar loops.
- The legacy job scripts use retired settings (`par-single`, the `ivybridge` constraint, `intel/17`
  modules, `/gws/nopw/j04` paths) and will not work.

## Storage (be careful)
| Where | Size / limit | Use |
|---|---|---|
| Home `~` | 100 GB quota (~9 GB used, 2026-11-09) | code, small outputs only; don't push it |
| `~/acrg/mattrigby` = `/gws/ssde/j25b/acrg/mattrigby` | ACRG GWS, **97% full (~700 GB free of 22 TB, 2026-11-09)** | larger outputs to keep |
| `/work/scratch-pw[4-7]/mattrigby` | max 100 TB in total; files **deleted 28 days after last access** | NAME working directories, restored met, intermediate files |

Met restore (official `MetRestore_JASMIN.sh`): Mk6 files are **copied and unzipped** into the met
folder (one day ≈ 3 GB gzipped, ≈ 12 GB unzipped, so ≈ 390 GB for a 30-day back-run). Mk10–12 are
linked, not copied. Use one shared met folder in scratch for all days of a month.

Run NAME in scratch, copy only what we keep (processed footprints, logs, inputs) to the GWS folder.
Avoid putting more than 100,000 files in one directory.

## Paths
| What | Path |
|---|---|
| NAME GWS (executables, UM met) | `~/name` = `/gws/ssde/j25a/name` |
| NAME 8.7.1 | `~/name/code/NAME-8.7.1` (`nameiii_64bit_par.exe`, `nameiii_eccodes_NetCDF_64bit_par.exe`) |
| UM topography | `/gws/ssde/j25a/name/data/UMTopogData` |
| ACRG GWS | `/gws/ssde/j25b/acrg` |

## UM global met coverage (archive, as of 2026-11-06)
| Version | Coverage |
|---|---|
| Mk10 (L59PT) | 2017-06-01 → 2022-05-04 12Z |
| Mk11 (L59PT) | 2022-04-01 → 2026-01-21 21Z |
| Mk12 (L59PT) | 2025-12-01 → ongoing (NAME 8.7.1 only) |

Met top ≈ 29 km (eta 0.364 × 80 km).

## Running NAME 8.7.1 (checked 2026-11-09, debug job 66266117)
Environment, as in the NAME team's demo scripts (`dev/scripts/prepare_debug_run.py` writes it):

```bash
module load jaspy
module load oneapi/compilers/24.2.0 oneapi/mpi/24.2.0
module load netcdf/intel2024.2.0/4.9.2 netcdf/intel2024.2.0/fortran/4.6.1
export OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_STACKSIZE=32m
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gws/ssde/j25a/name/lib/shumlib_2024.11.1/jasmin-lotus-ifx-gcc/lib
mpirun -np 1 /gws/ssde/j25a/name/code/NAME-8.7.1/Executables_Linux/nameiii_64bit_par.exe input.txt
```

Without the shumlib path NAME fails at start-up (`libshum_wgdos_packing.so: cannot open shared
object file`). Input-file lessons: NAME requires `Restart File Options` and `Multiple Case Options`
blocks; a line whose last item is blank must end with a comma; `Max # Field Reqs` must be well
above the number of requests (mixing ratio adds internal fields; we use 20,000 like legacy).
Output: `Fields_grid1_C1_T1_<YYYYMMDDHHMM>.txt` (read with `name_output.read_fields`) and
`Particles_C1.txt` (fate output, `name_output.read_name_fate_file`).

## Python
Use the project's uv environment (`uv run ...`). The legacy code's tests run under
`module load jaspy/3.12/v20250704` if you need to compare against it.
