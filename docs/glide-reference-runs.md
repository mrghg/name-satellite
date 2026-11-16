# NAME reference runs for GLIDE

NAME back-run footprints made to validate the GLIDE model, to GLIDE's specification
(`docs/reference-runs.md` in github.com/mrghg/glide; the version used is `spec.SPEC_VERSION`).
These runs are separate from the satellite pipeline. They reuse its NAME input writer
(`name_input.py`), job scripts (`lotus.py`) and met handling (`met.py`), but none of its sounding
or column code. Code: `src/name_satellite/glide_runs/`.

## What the spec asks for (summary; the spec is authoritative)

- 7 receptors (`spec.RECEPTORS`), one a column (Lamont), on two 0.25° domains, EUR and NAM
  (`spec.EUR`, `spec.NAM`); particles are removed when they leave the box or rise above 20 km agl.
- Hourly one-hour releases over two 14-day blocks in 2019 (`spec.BLOCKS`), 4,704 in all; every
  particle followed for 240 h.
- One footprint per release: mixing ratio in the 0-40 m layer, integrated over 240 h, converted
  to (mol/mol)/(mol m⁻² s⁻¹).
- Settings shared with FLEXPART take the same value (§7): sync time 300 s, deep convection on
  (`New`), boundary-layer depth 100-4,500 m, 6 particles/s (28/s for the column).
- Runs (§9): `main`; `seed1`-`seed5` (`Random` seed); 7 single-release runs with particle
  positions; `nomeso` (mesoscale motions off) and `group` (the group's site-template settings).

## Workflow

1. **Prepare a run** (`glide_runs.prepare`): one receptor, a list of hourly release times and a
   variant. Writes `plan.csv` (one NAME source per release, `glide_runs.plan`), `input.txt`,
   `job.sbatch`, `restore_met.sbatch` and `run.json`.
   ```
   uv run python -m name_satellite.glide_runs.prepare RUN_FOLDER MHD 2019-01-08T00 2019-01-08T23
   ```
2. **Run on LOTUS2**: submit `restore_met.sbatch`, then `job.sbatch` with
   `--dependency=afterok:<restore job>`. NAME is the ACRG-patched 8.7.1 build (its fate records
   carry each particle's index and why it stopped).
3. **Deliver** (to do: `glide_runs.deliver`): per receptor, block and run, the §8 footprint file
   and `_end` file, and the run record.

Raw output and deliverables go to `/gws/ssde/j25b/acrg/mattrigby/glide-reference-runs/`
(restricted: NAME output and the run record, which contains NAME's met declarations, never go
into a git repository). Tests go to `/work/scratch-pw5/$USER/name-satellite/glide-tests/`.

## The NAME input, part by part

- **Sources:** one per release, `Cuboid` in `m agl` with zero size (a point; dZ in metres), rate
  `6/sec`, from the end of the hour back to its start (backwards run).
- **Footprints:** one integrated `Mixing Ratio` request per source on the §2 grid, 0-40 m layer,
  one output time at the chronologically earliest end of the window covering every release
  (`name_input.run_timing`).
- **Particle fate records** (`Particles_C1.txt`): one per particle when it stops (left a side,
  left the top, or reached 240 h), for `n_particles`, `fraction_left_*` and the end records.
- **Receptor met** (`ReleaseMet_R1_C1.txt`): model ground height (`Topography`), surface pressure
  and boundary-layer depth at the receptor at the first release (§3, §8).
- **Particle positions** (single-release runs only, `Positions_C1.txt`): every particle each hour.

## Findings so far (details in `dev/notes/name-behaviour.md`)

- 16 Nov 2026, smoke test (job 2295092): all the spec's settings are accepted. A one-hour source
  at 6/s releases 21,601 particles. The 240 h limit is exact (a particle stops at that age).
  `Land Fraction` is not in the Mk10 met, so the receptor's surface type cannot be reported from
  NAME (§3 allows this: "if the model can report it").
- 16 Nov 2026, test 2 (job 2296900): a source in `Pa` spreads particles evenly in pressure, so
  each Lamont release is ONE source from the model surface pressure to 100 hPa (§5), not 28 layers.
- 16 Nov 2026, test 3 (job 2296902; §10 check 6): two identical sources in a `Fixed` run are not
  identical, so `main` uses `Fixed` (§9).
- 16 Nov 2026, memory probe (job 2301599): 5.2 M live particles (a full point-receptor block)
  need 5.1 GB; ~0.9 CPU-min per million particle-hours with convection on.
- 16 Nov 2026, §10 one-day check at Mace Head (jobs 2306225-8, folder
  `glide-reference-runs/check10-2026-11-16/` on the ACRG GWS, results in `check10_results.txt`,
  run record `run_record_NAME_check10.yaml`): checks 1, 2, 5 and 7 pass; 21,601 particles per
  release; two `Fixed` runs bit-identical; footprint totals median 3.9 m² s mol⁻¹ (below the
  group's earlier 5-17, not by a factor of ten); 60-100% of particles leave the EUR box through
  its sides (87% through the west edge), none through the top. Check 8: patched and stock builds
  give identical fate records but one footprint value differs in the 7th significant figure.
  Each 240 h one-day run took 21-31 min on 8 threads (about 5 GB).

## Still to do

- Lamont: hourly model surface pressure at the receptor (met-only run), then the column plan.
- Send the §10 files and run record to the GLIDE side; then the blocks (14 `main` runs), the
  repeat-seed, single-release, `nomeso` and `group` runs.
