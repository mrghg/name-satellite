# Implementation plan

Guiding principle: **first reproduce the legacy workflow and test against it; then fix bugs one
at a time as logged deviations; then add new features.** Agreed with Matt Rigby, 2026-11-06.

- **Deterministic steps** (release geometry, NAME input generation, reading output, processing)
  must reproduce legacy results exactly, or within floating-point rounding, from the same inputs.
- **NAME runs themselves** are stochastic. A re-run can only be compared statistically, using the
  spread between independent legacy runs of the same day as the noise yardstick (see
  [reference-data.md](reference-data.md)).

Status key: ☐ not started, ◐ in progress, ☑ done.

## Phase 0: Ground truth and repo ☑/◐
- ☑ Repo, uv, ruff, gitleaks, pytest, docs.
- ☑ Reference data for 2013-01-01 GOSAT-BRAZIL registered with checksums.
- ☑ Showed the legacy NAME sources can be regenerated exactly from the obs file plus the NAME surface
  pressure (`dev/scripts/reproduce_legacy_sources.py`).
- ☐ Rest of January 2013 raw output (from Matt's archive) for month-scale checks.
- ☐ Remaining open questions: see [decisions.md](decisions.md#open-questions).

## Phase 1: Satellite input (L2 → canonical sounding table)
- Reader for University of Leicester OCPR proxy XCH4: GOSAT-1 v10.0 and GOSAT-2 v2.0 (first
  supported input), plus ESA CCI v7.x for legacy validation only.
- Canonical sounding table with explicit, stable IDs (`exposure_id`-based), time, lat/lon, footprint
  size, pressure levels, pressure weights, averaging kernel, a priori, XCH4 and uncertainty,
  flags, product and version.
- Port the legacy filters and options: quality flag, bad pressure, land/glint mode, optional binning
  to the NAME grid. Test them against `acrg.satellite.gosat` on the same file.
- Write the observation netCDF in the legacy format (dims `time`, `lev`, `id`).

## Phase 2: Release planning
- Release layers from retrieval levels (legacy `midpoint_bounds`, `set_as_edges=False`). Golden
  test: reproduce the 425 sources in the legacy input file.
- **Pressure-based release ceiling.** Any ceiling below the NAME domain top is allowed. Validate
  that ceiling < domain Z max < met top (about 29 km).
- Bottom-layer anchoring: a decision is pending (legacy used a separate NAME surface-pressure
  run).
- Particle allocation per level (legacy: 30/s at or below `LevelSplit`, otherwise 3/s, over a 60 s
  release).
- An explicit "release plan" object (sounding ID, level, pressure bounds, particle rate, optional
  weight), so the future weighted-release option can slot in.

## Phase 3: NAME input generation
- Python templating of NAME input files (no ksh). Golden tests against the legacy
  `Input_files` (NAME 7.2 syntax), then target NAME 8.7.1.
- Met registry built from the actual archive coverage (Mk10 → Mk11 → Mk12) rather than hard-coded
  epoch boundaries. Choose the run length and met version from what exists.
- Domains as configuration (grid, dispersion domain, Z max).
- Time-integrated and hourly (HiTRes) output modes.
- Release-point met without the per-release grid explosion (the legacy method fails with
  "too many horizontal grids" for a normal day).

## Phase 4: Running on LOTUS2
- SLURM job generation (array per day), scratch handling, met restore, manifests, idempotent
  restarts, captured logs, and failure detection that parses NAME's logs.

## Phase 5: Reading NAME output
- Fields text files (time-integrated and hourly), particle fate files (legacy ACRG format and
  stock NAME 8.7.1 fate format), met output.
- Match by source/request name, never by column order. Check particle counts per source.

## Phase 6: Processing
- Port the satellite path of `NAME_processing.process()` as small pure functions:
  - unit conversion
  - footprint array
  - boundary-particle histograms
  - release met
  - column combination with pressure weights × averaging kernel and `max_level`
- Reproduce the benchmark `fp` file for 2013-01-01 exactly, then the full month.
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
