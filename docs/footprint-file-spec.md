# NAME satellite footprint file: draft specification

Draft, 2026-11-10 (Matt Rigby). Under review by the group: until it is agreed, the pipeline still
writes the legacy layout (`legacy_pipeline.write_legacy_fp_file`). Discussion copy:
https://claude.ai/artifact/DazSraoqEokkCNLXaPRmRr

## Purpose and scope

The new pipeline will write one netCDF footprint file per satellite product per month, keyed by
sounding ID, to replace the 2019 file `GOSAT-BRAZIL-column_SOUTHAMERICA_201301.nc`. It holds, for
every sounding, the column footprint, the boundary sensitivities and the met at the release point:
everything the inversion needs from NAME for that observation.

- **Readers:** the group's inversion code (today `acrg.name.footprints_data_merge` and the hbmcmc
  satellite path), and our own checks.
- **Writer:** `name-satellite`, after the NAME runs (workflow step 6).
- **Not in scope:** the observation file itself (already keyed by sounding ID, decision D18), hourly
  (HiTRes) footprints, and the `srr` format used by the newest inversions (a later, separate
  output).

## Design principles

Five rules, each fixing a problem found in the legacy file:

1. **Keyed by sounding ID.** The main dimension is `sounding_id` (the L2 exposure ID; for a bin,
   that of its earliest exposure). Time is an ordinary variable. GOSAT-1 v10 and GOSAT-2 times are
   whole minutes and repeat (X22), so time cannot identify a sounding (D18).
2. **One satellite product per file.** GOSAT-1 v10.0 and GOSAT-2 v2.0 are never mixed (D17); the
   product and version are in the file name and attributes.
3. **No lossy rounding.** Values are stored as float32 without `least_significant_digit`. The
   legacy rounding lowered column totals by 0.4% on average, up to 1.4% (X7).
4. **Units on every variable**, stated correctly (legacy labelled °C as K, X17, and Pa as m, X8).
5. **Provenance.** Each file records the code version, NAME version, met, run settings and the list
   of deviations from legacy applied, so a file can be traced and reproduced.

## File layout

One file per satellite product, domain and month, as legacy, named for example
`GOSAT1-v10.0_BRAZIL-column_SOUTHAMERICA_201301.nc` (legacy:
`GOSAT-BRAZIL-column_SOUTHAMERICA_201301.nc`).

| Dimension | Size (Jan 2013 example) | Coordinate and units |
| --- | --- | --- |
| `sounding_id` | 678 | sounding ID (text), sorted by time, then ID |
| `lat` | 357 | cell-centre latitude, degrees_north |
| `lon` | 190 | cell-centre longitude, degrees_east |
| `pressure_bin` | ~20 (to decide) | bin centre, Pa, with `pressure_bin_bounds` (D14) |
| `release_level` | 16 | 1 = surface layer (16 modelled levels for the GOSAT products, D21) |

Large arrays keep legacy's order, sounding last (e.g. `fp(lat, lon, sounding_id)`), so the
inversion's reshapes still work. They are zlib-compressed, chunked one sounding per chunk. Lossless
storage is 2.6 times larger (1 Jan 2013: 2.14 MB against 0.81 MB), so a month is about 65 MB
against legacy's 25 MB.

## Variables

Every column quantity is already weighted by pressure weight × averaging kernel and summed over the
levels NAME models, as now (Ganesan et al. 2017, eq. 5).

| Variable | Dimensions | Units | Meaning |
| --- | --- | --- | --- |
| `time` | sounding_id | datetime (UTC) | measurement time from the L2 file |
| `fp` | lat, lon, sounding_id | (mol/mol) / (mol m⁻² s⁻¹) | column footprint |
| `bs_n`, `bs_s` | pressure_bin, lon, sounding_id | 1 | column-weighted share of particles leaving through the north / south face |
| `bs_e`, `bs_w` | pressure_bin, lat, sounding_id | 1 | the same for the east / west faces |
| `bs_top` | lat, lon, sounding_id | 1 | leaving through the domain top (D13) |
| `fraction_inside` | sounding_id | 1 | still inside the domain at the maximum travel time |
| `column_weight` | release_level, sounding_id | 1 | pressure weight × averaging kernel for each modelled level |
| `release_pressure` | release_level, sounding_id | Pa | centre of each release layer |
| `release_thickness` | release_level, sounding_id | Pa | thickness of each release layer |
| `release_lat`, `release_lon` | sounding_id | degrees_north / degrees_east | centre of the release box |
| `release_dlat`, `release_dlon` | sounding_id | degrees | full size of the release box (widened for bins) |
| `n_soundings` | sounding_id | 1 | L2 soundings averaged into the entry (1 if not binned) |
| `temperature` | sounding_id | K | air temperature at the lowest release |
| `pressure` | sounding_id | Pa | air pressure at the lowest release |
| `surface_pressure` | sounding_id | Pa | model surface pressure at the release point |
| `wind_speed` | sounding_id | m s⁻¹ | wind speed at the lowest release |
| `wind_direction` | sounding_id | degrees | direction the wind blows from |
| `boundary_layer_depth` | sounding_id | m | at the release point |

Check sum per sounding: the boundary sensitivities, `bs_top` and `fraction_inside` add up to
Σ `column_weight`, so every particle is counted exactly once. That rule depends on the NAME team's
answer about duplicate fate records (X20).

## Global attributes and provenance

Enough to trace any file back to its inputs and code, and to check it matches its observation file:

| Attribute | Example |
| --- | --- |
| `product`, `product_version` | UoL GOSAT OCPR, 10.0 |
| `l2_files` | the L2 file names used, one per day |
| `domain`, `site` | SOUTHAMERICA, GOSAT-BRAZIL |
| `max_level` | 16 (levels modelled by NAME, D21; the inversion checks it against the observation file) |
| `name_version` | NAME III 8.7.1 |
| `met` | UMG_Mk6PT |
| `max_travel_hours`, `integration_window` | 720, all_releases (X18) |
| `particles_per_release` | rule, e.g. 30/s for levels 1-12, 3/s above, 60 s |
| `random_seed` | Fixed |
| `deviations_applied` | e.g. X7, X18, X21 (IDs from docs/deviations.md) |
| `software`, `software_commit` | name-satellite, git commit |
| `created`, `Conventions` | ISO time, CF-1.8 |

## Differences from the legacy file

The footprint keeps its name and units; most other changes are renames, unit fixes and the new
boundary sensitivities.

| Legacy (2019) | New | Change |
| --- | --- | --- |
| dimension `time` | `sounding_id`, with `time` as a variable | key is the sounding ID (D18) |
| `fp` | `fp` | no lossy rounding (X7); window covers every sounding (X18) |
| `particle_locations_n/e/s/w` (height, along edge, time) | `bs_n/e/s/w` (pressure_bin, along edge, sounding_id) | every particle once (D13, X6, X14); pressure bins, not m (D14) |
| none | `bs_top`, `fraction_inside` | particles leaving through the top, or still inside (D13) |
| `height` (m, no reference) | `pressure_bin` (Pa) with bounds | D14 |
| `temperature` (°C labelled K) | `temperature` (K) | X17 |
| `pressure` (hPa) | `pressure` (Pa); plus `surface_pressure` | SI units throughout |
| `PBLH` | `boundary_layer_depth` (m) | rename |
| `release_lat`, `release_lon` | same, plus `release_dlat`, `release_dlon` | release box stored |
| none | `column_weight`, `release_pressure`, `release_thickness`, `n_soundings` | how each column was built |
| `lev` = "column" | removed | always a column |
| attributes `max_level`, `author`, `created` | full provenance (previous section) | |

## Changes needed in the inversion code

The one change that is not optional: observations and footprints must be joined on `sounding_id`.
Today `footprints_data_merge` joins them by time, forward-filling within 60 s (`acrg/name/name.py`
~1013-1037, `combine_datasets` ~723). With GOSAT-1 v10 and GOSAT-2 times repeating to the minute
(X22), that can silently attach the wrong footprint.

1. Join observations and footprints exactly on `sounding_id` (both files carry it), and fail on any
   unmatched ID.
2. Use `sounding_id` wherever `time` is the observation dimension (e.g. the `H` reshape ~1305, the
   dimension order ~1039). Keep `time` for selecting fluxes by month.
3. Boundary conditions (~1208, `bc_sensitivity` ~1400): use `bs_n/e/s/w` on `pressure_bin`, with
   the global-model curtains interpolated to pressure (D14), and decide how `bs_top` and
   `fraction_inside` enter (D13).
4. Read `temperature` in K and `pressure` in Pa (no longer °C and hPa).
5. Keep the `max_level` check between observation and footprint files.

For a transition period the writer could also add legacy-named variables, so the current code runs
unchanged on the new files (open question below).

## Open questions for the group

- [ ] Monthly files as now, or one file per day (simpler to rerun a single day)?
- [ ] Pressure bins for the boundary sensitivities: how many, and their edges (D14)?
- [ ] `bs_top` on the full lat × lon grid is as large as `fp`; is a coarser grid, or a few latitude
  bands, enough for the stratospheric boundary term?
- [ ] Write legacy-named copies (`particle_locations_*`, `PBLH`, hPa pressure) during a transition,
  or switch the inversion code in one go?
- [ ] Store per-level footprints as well as the column (16 times larger), e.g. to test other
  averaging kernels later?
- [ ] Who updates the inversion code, and against which branch?
- [ ] Soundings sharing NAME releases (D20): add `release_group` (sounding_id, text: the group's
  ID, the sounding ID of its earliest member) so the inversion can treat a group's soundings
  together (e.g. correlated errors, or its own averaging)? Their footprints are not independent
  samples: they come from the same particles.
