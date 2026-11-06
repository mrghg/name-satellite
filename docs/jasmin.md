# Running on JASMIN

## Where to run what
- **Sci / login node:** editing, `uv sync`, unit tests, single-day processing, inspecting output.
  Keep runs short and modest in memory; these nodes are shared.
- **LOTUS2 (SLURM):** all NAME runs and anything heavy (month-scale processing, many days).
  Submit with `sbatch`; never run NAME interactively.

## LOTUS2 notes
- Partitions include `standard`, `debug`, `special`, `interactive`. The legacy job scripts use
  retired settings (`par-single`, the `ivybridge` constraint, `intel/17` modules, `/gws/nopw/j04`
  paths) and will not work.
- Account / QoS for this project: **TBD** (open question 4 in decisions.md).
- Scratch: `/work/scratch-pw[4-7]/$USER` (fast, not backed up, purged). Put NAME working
  directories and restored met here.

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
