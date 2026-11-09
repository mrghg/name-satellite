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

## Python
Use the project's uv environment (`uv run ...`). The legacy code's tests run under
`module load jaspy/3.12/v20250704` if you need to compare against it.
