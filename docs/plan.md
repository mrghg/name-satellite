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

## Phase 0: Ground truth and repo ☑/◐
- ☑ Repo, uv, ruff, gitleaks, pytest, docs.
- ☑ Reference data for 2013-01-01 GOSAT-BRAZIL registered with checksums.
- ☑ Showed the legacy NAME sources can be regenerated exactly from the obs file plus the NAME surface
  pressure (`dev/scripts/reproduce_legacy_sources.py`).
- ☑ January 2013 raw output of both legacy runs (2018 conc, 2019 MR). The legacy processing port
  reproduces the production file for all 31 days: 26 bit-identical, 5 with one footprint cell one
  rounding step off ([reference-data.md](reference-data.md#new-code-vs-legacy-2026-11-07)).
- ☐ Remaining open questions: see [decisions.md](decisions.md#open-questions).

## Phase 1: Satellite input (L2 → canonical sounding table) ◐
- ☑ `observations.py` (2026-11-09): reader for UoL OCPR GOSAT-1 v10.0, GOSAT-2 v2.0 and ESA CCI
  v7.x; filters (quality, undefined pressure, level order, land/glint, lat/lon); legacy binning;
  legacy observation-file writer. Identical to `acrg.satellite.gosat` on v7.0 files for 1 and 7 Jan
  2013 (the legacy code cannot read v10/GOSAT-2). Findings: X22 (minute times), X23 (binning
  quirks), X24 (bad-pressure filter).
- ☑ Soundings keyed by sounding ID (D18): `observations.to_obs`, `sounding_id` in release plans.
- ◐ Footprint files keyed by `sounding_id` too (improved processing): draft format in
  [footprint-file-spec.md](footprint-file-spec.md), under review by the group; the inversion's
  observation reading to be switched to it (with the group).
- ☐ Reference data for a durable test against legacy (v7.0 L2 files + legacy output, to `data/`).
- Reader for University of Leicester OCPR proxy XCH4: GOSAT-1 v10.0 and GOSAT-2 v2.0 (first
  supported input), plus ESA CCI v7.x for legacy validation only.
- Canonical sounding table with explicit, stable IDs (`exposure_id`-based), time, lat/lon, footprint
  size, pressure levels, pressure weights, averaging kernel, a priori, XCH4 and uncertainty,
  flags, product and version.
- Port the legacy filters and options: quality flag, bad pressure, land/glint mode, optional binning
  to the NAME grid. Test them against `acrg.satellite.gosat` on the same file.
- Write the observation netCDF in the legacy format (dims `time`, `lev`, `id`).

## Phase 2: Release planning ◐
- ☑ Release layers from retrieval levels (legacy `midpoint_bounds`, `set_as_edges=False`):
  `releases.release_layers_pa`. Golden tests reproduce all legacy sources (1 Jan MR run, and all 31
  days / 678 soundings / 11,526 sources of the conc run): levels 3–17 from the obs file alone;
  levels 1–2 also once the legacy NAME surface pressure (inferred) is supplied.
- ☑ Legacy particle-rate rule: `releases.particles_per_second`.
- **Pressure-based release ceiling.** Any ceiling below the NAME domain top is allowed. Validate
  that ceiling < domain Z max < met top (about 29 km).
- Bottom-layer anchoring: a decision is pending (legacy used a separate NAME surface-pressure
  run).
- Particle allocation per level (legacy: 30/s at or below `LevelSplit`, otherwise 3/s, over a 60 s
  release).
- An explicit "release plan" object (sounding ID, level, pressure bounds, particle rate, optional
  weight), so the future weighted-release option can slot in.

## Phase 3: NAME input generation ◐
- ☑ Python writer for NAME 8.7.1 input files (`name_input.py`); sources round-trip exactly and
  match the legacy sources; legacy timing reproducible (`integration="legacy_2019"`).
- ☑ NAME 8.7.1 runs our input file (LOTUS2 debug job 66266117: 2 releases, 3 h, Mk6 met restored
  automatically; fields and fate files read by our readers).
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
- ◐ Stage 2 validation (1–7 Jan 2013; `dev/scripts/prepare_stage2.py`, `process_stage2.py`):
  daily runs for 2–7 Jan with the legacy window and legacy release boxes, compared with
  production (7 Jan: 1,088 sources), and all 7 days in ONE run (2,720 sources, corrected
  window) to test batching (`batching.py`). Jobs 66323615 (met restore) and 66323616–22.
- ☑ Met registry (`met.py`): official 8.7.1 version dates, file names from NAME's met
  definitions, archive check. ☐ Runs spanning a changeover (currently refused).
- Domains as configuration (grid, dispersion domain, Z max).
- Time-integrated and hourly (HiTRes) output modes.
- ☑ Release-point met from the main run (`release_met.py`; one location set and pressure array
  per ≤86 soundings, one time grid and output file per sounding; up to 199 soundings per run).
  Checked against the legacy met run on 1–2 Jan 2013: footprints and histograms bit-identical,
  met within NAME's print precision; underground releases get the model surface pressure, as in
  legacy (`dev/scripts/prepare_release_met_test.py`, `compare_release_met.py`).

## Phase 4: Running on LOTUS2
- SLURM job generation (array per day), scratch handling, met restore, manifests, idempotent
  restarts, captured logs, and failure detection that parses NAME's logs.

## Phase 5: Reading NAME output ◐
- ☑ Time-integrated fields text files (`name_output.read_fields`), legacy ACRG particle files
  (`read_acrg_particle_file`) and legacy release-point met (`legacy.read_legacy_met_folder`), checked
  against the source list for all 31 days of January 2013.
- ☑ Stock NAME 8.7.1 fate format (`read_name_fate_file`).
- ☐ Hourly fields; new release-met approach.
- Match by source/request name, never by column order. Check particle counts per source.

## Phase 6: Processing ◐ (legacy reproduction ☑)
- Port the satellite path of `NAME_processing.process()` as small pure functions:
  - unit conversion
  - footprint array
  - boundary-particle histograms
  - release met
  - column combination with pressure weights × averaging kernel and `max_level`
- ☑ Reproduce the production `fp` file for 2013-01-01 exactly: every variable bit-identical
  (`legacy_pipeline.process_day_legacy` + `write_legacy_fp_file`). The full month can't be checked
  exactly because only 1 Jan of the mixing-ratio raw output survives.
- Output both the legacy `fp` format (used by the current satellite inversion code) and `srr` (used
  by the newest inversions), as decided.

## Phase 7: Fix legacy bugs (each a logged deviation)
- Candidates are listed in [deviations.md](deviations.md). Each fix comes with a test that
  quantifies its effect on the 2013-01-01 and January 2013 footprints.

## Phase 8: New capability
- NAME 8.7.1 with Mk11/Mk12 met; end-to-end run of one recent day on LOTUS2.
- Hourly mode validated against time-integrated (sum over hours = integrated).
- Other instruments (GOSAT-2 first, then TROPOMI and OCO-2) via the canonical table.

## Phase 9: Architecture evaluation (after the baseline is validated)
- (a) Release particles weighted by pressure weight × averaging kernel, using Source Groups so
  there is one column field per sounding. Validate against Σ wᵢ·fpᵢ from per-level runs.
- (b) NAME NetCDF field output (Output Route `N`, NetCDF executable) instead of text. Measure
  size and speed on a real day.
