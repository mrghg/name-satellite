# Implementation plan

Guiding principle: **first reproduce the legacy workflow and test against it; then fix bugs one
at a time as logged deviations; then add new features.** Agreed with Matt Rigby, 2026-11-06.

- **Deterministic steps** (release geometry, NAME input generation, reading output, processing)
  must reproduce legacy results exactly, or within floating-point rounding, from the same inputs.
- **NAME runs themselves** are stochastic. A re-run can only be compared statistically, using the
  spread between independent legacy runs of the same day as the noise yardstick (see
  [reference-data.md](reference-data.md)).

Code is written for **clarity rather than perfection** (readable by scientists; workflow described
in docstrings and [workflow.md](workflow.md)).

Status key: ☐ not started, ◐ in progress, ☑ done.

## Phase 0: Ground truth and repo ☑
- ☑ Repo, uv, ruff, gitleaks, pytest, docs.
- ☑ Reference data for 2013-01-01 GOSAT-BRAZIL registered with checksums.
- ☑ Showed the legacy NAME sources can be regenerated exactly from the obs file plus the NAME surface
  pressure (`dev/scripts/reproduce_legacy_sources.py`).
- ☑ January 2013 raw output of both legacy runs (2018 conc, 2019 MR). The legacy processing port
  reproduces the production file for all 31 days: 26 bit-identical, 5 with one footprint cell one
  rounding step off ([reference-data.md](reference-data.md#new-code-vs-legacy-2026-11-07)).
- Open questions: see [decisions.md](decisions.md#open-questions).

## Phase 1: Satellite input (L2 → canonical sounding table) ☑
- ☑ `observations.py` (2026-11-09): reader for UoL OCPR GOSAT-1 v10.0, GOSAT-2 v2.0 and ESA CCI
  v7.x; sounding table with explicit IDs, product and version; filters (quality, undefined
  pressure, level order, land/glint, lat/lon); legacy binning; legacy observation-file writer.
  Identical to `acrg.satellite.gosat` on v7.0 files for 1 and 7 Jan 2013 (reference data
  `gosat_cci_v7.0_201301_legacy_obs`, `tests/test_observations.py`; the legacy code cannot read
  v10/GOSAT-2). Findings: X22 (minute times), X23 (binning quirks), X24 (bad-pressure filter).
- ☑ Soundings keyed by sounding ID (D18): `observations.to_obs` / `write_obs_file`, `sounding_id`
  in release plans.
- ☑ L2 files downloaded from CEDA on demand (`l2_files.py`).

## Phase 2: Release planning ◐
- ☑ Release layers from retrieval levels (legacy `midpoint_bounds`, `set_as_edges=False`):
  `releases.release_layers_pa`. Golden tests reproduce all legacy sources (1 Jan MR run, and all 31
  days / 678 soundings / 11,526 sources of the conc run): levels 3–17 from the obs file alone;
  levels 1–2 also once the legacy NAME surface pressure (inferred) is supplied.
- ☑ Legacy particle-rate rule (30/s up to `LevelSplit` = 12, 3/s above, 60 s):
  `releases.particles_per_second`.
- ☑ Release plan table (`releases.build_release_plan`): one row per release with sounding ID,
  level, position and box, pressure centre and thickness, times and particle rate. ☐ Optional
  per-release weight (deferred option (a)).
- ☑ Bottom layer (D19): new runs scale every level to the model surface pressure, from a met-only
  NAME run before the main run (`surface_pressure.py`, chained by `prepare.py`); validation runs
  use the legacy NAME surface pressure (D16). ☑ Effect against the retrieval's own levels on
  1 Jan 2013 (30-day runs, same seed, jobs 66679485/6/8, `dev/scripts/prepare_d19_comparison.py`,
  `compare_d19.py`; 46 soundings): column footprint total median +2.6% (5-95% −2.9 to +10.5%),
  within 2° of the sounding +6.7% (−1.7 to +32%); shape differences at the particle-noise level
  (16x coarsened: median 5%). 1-h test (jobs 66678616/8): the
  main run's surface pressure equals the pre-run's exactly; no level-1 centre below ground.
- ☑ Release groups (D20, `release_groups.py`): soundings in the same output-grid cell and within
  60 s share one set of releases; each sounding keeps its own values (`release_groups.csv`).
  ☐ Processing: column footprint per sounding from its group's per-level footprints, retrieval
  level j ↔ release level j+1 by number (not legacy nearest pressure).
- ☑ Modelled levels (D21, X10): levels wholly inside a 20 km agl domain (16 for the GOSAT
  products), prior above; per-sounding checks of the domain top (from the met-only run) and the
  modelled column weight (≥ 0.85). ☐ Measure the effect against 17 levels / 19 km.
- ☑ Particle rates (D22, `releases.PARTICLE_RATES`): 30/s for levels 1–13, 10/s for level 14,
  3/s for 15–16; from the stage 1 seed pair (`dev/scripts/analyse_level_noise.py`).

## Phase 3: NAME input generation and validation ☑ (open items below)
- ☑ Python writer for NAME 8.7.1 input files (`name_input.py`); sources round-trip exactly and
  match the legacy sources; legacy timing reproducible (`integration="legacy_2019"`).
- ☑ NAME 8.7.1 runs our input file (LOTUS2 debug job 66266117).
- ☑ Stage 1 validation (1 Jan 2013, all 425 releases, 30 days, Mk6; LOTUS2 jobs 66271462/3, 86 min
  each on 4 threads). Noise test: run B2 is run B with a different random seed (job 66320292).
  Run B (legacy window) differs from production by no more than B differs from B2, on every
  measure: column totals (scatter 4.1% vs 4.3%), footprint shapes at every coarsening, edge
  histograms, and per-release totals at every level. Column totals: B/production +0.5 ± 0.8%,
  B2/production −0.4 ± 0.8%. Open, minor: the upper (3 particles/s) releases total 10–14% more
  than production in both B and B2; these releases make up a median 4% of a column's footprint, so
  ≤0.5% on columns, within noise, and one production run can't separate this from its own noise.
  Run A confirms X18. Script: `dev/scripts/process_stage1.py`; see [deviations.md](deviations.md)
  X18, X20.
- ☑ Stage 2 validation (2–7 Jan 2013, one run per day, legacy window and boxes;
  `dev/scripts/prepare_stage2.py`, `process_stage2.py`): column totals within noise every day
  (median ratio 0.99–1.01; 160 soundings with stage 1: scatter 3.4% vs the 4.3% seed yardstick);
  7 Jan shapes differ a little more than the 1 Jan yardstick, judged minor (no second seed run).
- ☑ Batching (`batching.py`): several days in one run, split afterwards. The first 7-day run failed
  at the end on NAME's output line limit (fixed: field requests in groups of ≤ 2,000); the fix was
  checked with a 1-hour-travel run of all 2,720 releases. The full 8-h rerun was skipped: batching
  saves little compute (particle work dominates), so daily runs are the default.
- ☑ Met registry (`met.py`): official 8.7.1 version dates, file names from NAME's met
  definitions, archive check. ☐ Runs spanning a met changeover (currently refused; X2).
- ☑ Domains as configuration (`domains.py`: SOUTHAMERICA, GOSAT-BRAZIL). Domains are being
  re-evaluated (open question 11).
- ☑ Release-point met from the main run (`release_met.py`; one location set and pressure array
  per ≤86 soundings, one time grid and output file per sounding; up to 199 soundings per run).
  Checked against the legacy met run on 1–2 Jan 2013: footprints and histograms bit-identical,
  met within NAME's print precision (X21).
- ☑ Fate output with met columns (`RunSettings.fate_met`): height in m agl plus pressure.
- ☐ Hourly (HiTRes) output mode (D1).

## Phase 4: Running on LOTUS2 ◐
- ☑ `prepare.py` (2026-11-10): a month of runs from current L2 data (download, select, observation
  file, release plan, NAME input, jobs, met restore, `submit.sh`). Tested: January 2013 GOSAT-1 v10
  prepared (30 days with data, 1,285 soundings, 21,845 releases: 1.9 times legacy, because of no
  binning and more good v10 soundings); 1 Jan run with 1-h travel read back and combined into 46
  column footprints. No full month run yet: domains first (open question 11).
- ☑ Job scripts: one per day, met restore as its own job, `afterok` dependencies; a job fails if
  NAME exits non-zero or writes to its error file.
- ☐ Array jobs, run manifests, restartable (idempotent) months, collected logs.
- ☑ Patched NAME 8.7.1 build with the particle index and fate code in the fate output (X20;
  `met.PATCHED_NAME_DIR`, outside the repository). X20 reported to the NAME team (2026-11-11).

## Phase 5: Reading NAME output ☑ (hourly to do)
- ☑ Time-integrated fields text files (`name_output.read_fields`; several output groups:
  `read_fields_folder`), keyed by source name, with grid checks.
- ☑ Legacy ACRG particle files (`read_acrg_particle_file`) and legacy release-point met
  (`legacy.read_legacy_met_folder`), checked against the source list for all 31 days of Jan 2013.
- ☑ Stock NAME 8.7.1 fate files (`read_name_fate_file`), with or without met columns.
- ☑ Release-point met from the main run (`release_met.read_release_met`).
- ☐ Hourly fields.

## Phase 6: Processing ◐ (legacy reproduction ☑)
- ☑ Satellite path of `NAME_processing.process()` ported as small functions (`footprints.py`,
  `boundary_histograms.py`, `column.py`, `legacy_pipeline.py`): reproduces the production `fp` file
  for all of January 2013 (26 days bit-identical, 5 with one cell one rounding step off).
- ☑ Lossless storage option (X7).
- ◐ New footprint file format keyed by `sounding_id`: draft in
  [footprint-file-spec.md](footprint-file-spec.md), under review by the group (supersedes D4 once
  agreed). Then: processing of new runs into it, and the inversion's reading switched to
  `sounding_id` (with the group).
- ☐ Improved boundary sensitivities (D13, D14), after X20.
- ☐ `srr` output (D4).

## Phase 7: Fix legacy bugs (each a logged deviation)
- Candidates are listed in [deviations.md](deviations.md). Each fix comes with a test that
  quantifies its effect on the 2013-01-01 and January 2013 footprints.

## Phase 8: New capability
- NAME 8.7.1 with Mk11/Mk12 met; end-to-end run of one recent day on LOTUS2.
- Hourly mode validated against time-integrated (sum over hours = integrated).
- Other instruments via the canonical table: ☑ GOSAT-2 reader; ☐ TROPOMI, OCO-2.

## Phase 9: Architecture evaluation (after the baseline is validated)
- (a) Release particles weighted by pressure weight × averaging kernel, using Source Groups so
  there is one column field per sounding. Validate against Σ wᵢ·fpᵢ from per-level runs.
- (b) NAME NetCDF field output (Output Route `N`, NetCDF executable) instead of text. Measure
  size and speed on a real day.
