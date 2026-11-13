# NAME satellite footprint file: draft specification

Draft, revised 2026-11-12 (Matt Rigby). Under review by the group: until it is agreed, the pipeline
still writes the legacy layout (`legacy_pipeline.write_legacy_fp_file`). Discussion copy:
https://claude.ai/artifact/DazSraoqEokkCNLXaPRmRr

## Purpose and scope

The new pipeline will write one netCDF footprint file per satellite product per month. It holds, for
every sounding, the column footprint, the boundary sensitivities and the met at the release point:
everything the inversion needs from NAME for that observation. It replaces the 2019 file
`GOSAT-BRAZIL-column_SOUTHAMERICA_201301.nc`.

The layout follows the group's current **site** footprint files (written by today's
NAME_processing; examples `MHD-10magl_NAME_UMG_EUROPE_inert_202501.nc` and, with time-resolved
footprints, `MHD-10magl_NAME_UMG_EUROPE_co2_202301.nc`), so that the inversion code can read
satellite and site footprints in the same way. It differs from them only where a satellite column
needs it (decision D23).

- **Readers:** the group's inversion code, and our own checks.
- **Writer:** `name-satellite`, after the NAME runs (workflow step 6).
- **Not in scope yet:** time-resolved footprints (`srr_time_resolved`, `srr_residual`, dimension
  `resolution`, as in the site files). They need hourly NAME output and will be specified once the
  rest of the file is agreed. The observation file itself is already keyed by sounding ID (D18).

## Design principles

1. **Same names and layout as the site files** wherever the meaning is the same: `srr`,
   `latitude`, `longitude`, `height`, `particle_locations_*`, `mean_age_particles_*`, the CF met
   names, the observation dimension first.
2. **Keyed by sounding ID, not time.** GOSAT-1 v10 and GOSAT-2 times are whole minutes and repeat
   (X22), so time cannot identify a sounding (D18). The observation dimension is `sounding_id`;
   `time` is an ordinary variable.
3. **Footprints are looked up through an index**, so soundings can later share a footprint without
   a format change. At first every sounding has its own (next section).
4. **One satellite product per file.** GOSAT-1 v10.0 and GOSAT-2 v2.0 are never mixed (D17).
5. **No biased rounding.** Legacy rounding to fixed decimal places lowered column totals by 0.4% on
   average, up to 1.4% (X7). The large arrays are rounded to a fixed number of *significant bits*
   instead, which never zeroes small values and leaves totals unbiased (Storage section).
6. **Units on every variable**, and **provenance**: each file records the code version, NAME
   version, met, run settings and the deviations from legacy applied.

## Soundings and footprints

Soundings in the same grid cell and minute share NAME releases (D20), but each sounding still gets
its own column footprint, built from its group's per-level footprints with its own pressure weights
and averaging kernel. The file separates the two:

```
dimensions: sounding_id = 1285 ; footprint = 1285 ; ...
int   footprint_index(sounding_id) ;      // which footprint each sounding uses (0-based)
float srr(footprint, latitude, longitude) ;
```

A reader gets each sounding's footprint with `ds.srr.isel(footprint=ds.footprint_index)`.

- **`footprint` = `sounding_id`** (`footprint_index` = 0, 1, 2, ...): every footprint is exact for
  its sounding. **Decided (D26, 2026-11-13): soundings of a release group do not share a
  footprint.** The measurement behind that decision follows.
- Considered: the soundings of a release group could share one footprint, after
  we have measured the error that introduces (members differ only through their weights; compare
  cell by cell, in totals, and in modelled XCH4 with a prior flux field, against the particle
  noise). The shared footprint would be labelled by the group's earliest sounding. Sharing would
  roughly halve `srr` (January 2013: 1,285 soundings in 660 groups).
  **Measured on 1 Jan 2013** (GOSAT-1 v10, 30-day run, 36 soundings in 14 multi-member groups;
  `dev/scripts/method_studies/measure_grouping_error.py`), sharing the group-mean-weight footprint (= the mean
  of the members' own) against each member's own:
  - `srr` total: median 0.21%, 90% 0.63%, max 0.90%; cell by cell (L1): median 0.41%, max 1.1%.
    About 15 times smaller than the particle noise in column totals (seed pair: ~3% per run), and
    zero on average within a group. With an illustrative uniform flux of 1e-8 mol/m²/s (own
    enhancements 4–55 ppb): median 0.02 ppb, max 0.42 ppb.
  - Boundary: the shared footprint's particle fractions add up to the group's column weight, not
    the sounding's (difference up to 2.4e-4): with a uniform 1850 ppb background, median 0.13 ppb,
    max 0.45 ppb. Rescaling the shared footprint by Σw_s / Σw̄ removes this exactly and leaves the
    other numbers unchanged (boundary arrays L1 median 0.18%).
  - Saving: about half of `srr` and the boundary arrays (~75 MB → ~40 MB for a GOSAT-1 month).
- The same index would serve legacy-style binning, or per-level footprints, if ever wanted.
- `release_group(sounding_id)` records which soundings share NAME releases (the group ID is the
  sounding ID of its earliest member). Their footprints come from the same particles, so they are
  not independent samples.

## File layout

One file per satellite product, domain and month, named like the site files
(`{site}_{model}_{met}_{domain}_{species}_{YYYYMM}.nc`), for example
`GOSAT1v10-BRAZIL-column_NAME_UMG_SOUTHAMERICA_inert_201301.nc` (exact form to agree).

| Dimension | Size (Jan 2013 example) | Coordinate and units |
| --- | --- | --- |
| `sounding_id` | 1,285 | sounding ID (text), sorted by time, then ID |
| `footprint` | 1,285 (= `sounding_id` at first) | `footprint_sounding_id`: the sounding whose weights built it |
| `latitude` | 357 | cell-centre latitude, degrees_north (double, as the site files) |
| `longitude` | 190 | cell-centre longitude, degrees_east |
| `height` | 20 | bin midpoint, m above ground, 1 km bins as the site files (D24); `height_bounds` |
| `release_level` | 16 | 1 = surface layer (16 modelled levels for the GOSAT products, D21) |

Large arrays put the observation dimension first, as the site files do (e.g.
`srr(footprint, latitude, longitude)`); the legacy satellite file put it last.

## Variables

Every column quantity is already weighted by pressure weight × averaging kernel and summed over the
levels NAME models, as now (Ganesan et al. 2017, eq. 5).

**On `footprint`** (column-weighted, so specific to the sounding whose weights built them):

| Variable | Dimensions | Units | Meaning | Site file |
| --- | --- | --- | --- | --- |
| `srr` | footprint, latitude, longitude | (mol/mol)/(mol/m2/s) | column footprint; `loss_lifetime_hrs` = -9 (inert) | same |
| `particle_locations_n`, `_s` | footprint, height, longitude | 1 | column-weighted share of particles leaving through the north / south face | same |
| `particle_locations_e`, `_w` | footprint, height, latitude | 1 | the same for the east / west faces | same |
| `mean_age_particles_n/e/s/w` | as above | hrs | mean travel time of those particles | same |
| `particle_locations_top` | footprint, latitude, longitude | 1 | leaving through the domain top (D13) | new |
| `mean_age_particles_top` | footprint, latitude, longitude | hrs | their mean travel time | new |
| `particle_fraction_inside` | footprint | 1 | still inside the domain at the maximum travel time (D13) | new |
| `footprint_sounding_id` | footprint | text | sounding whose weights built the footprint | new |

**On `sounding_id`:**

| Variable | Dimensions | Units | Meaning | Site file |
| --- | --- | --- | --- | --- |
| `time` | sounding_id | seconds since month start | measurement time from the L2 file (whole minutes, rounded down) | the dimension `time` |
| `footprint_index` | sounding_id | 1 | index into `footprint` | new |
| `release_group` | sounding_id | text | sounding ID of the group's earliest member (D20) | new |
| `air_temperature` | sounding_id | K | at the lowest release | same |
| `air_pressure` | sounding_id | hPa | at the lowest release | same |
| `wind_speed` | sounding_id | m s-1 | at the lowest release | same |
| `wind_from_direction` | sounding_id | degree | at the lowest release | same |
| `atmosphere_boundary_layer_thickness` | sounding_id | m | at the release point | same |
| `release_lat`, `release_lon` | sounding_id | degrees_north / degrees_east | centre of the group's release box | same |
| `release_dlat`, `release_dlon` | sounding_id | degrees | full size of the release box | new |
| `surface_pressure` | sounding_id | hPa | model surface pressure at the release point (D19) | new |
| `column_weight` | sounding_id, release_level | 1 | pressure weight × averaging kernel for each modelled level | new |
| `release_pressure` | sounding_id, release_level | hPa | centre of each release layer | replaces `release_height` |
| `release_thickness` | sounding_id, release_level | hPa | thickness of each release layer | new |
| `particles_released` | sounding_id, release_level | 1 | NAME particles released at each level, counted from the fate records | new |

The met is that of the group's shared release, so it is the same for all members of a group.

**Check sum per footprint:** the four `particle_locations_*`, `particle_locations_top` and
`particle_fraction_inside` add up to Σ `column_weight` of the sounding that built it, so every
particle is counted exactly once (D13, using the patched NAME's particle index, X20). In the site
files the four faces sum to 1; in the legacy satellite file they sum to the column weight (median
0.960), as here.

## Boundary heights

`height` is the height of the particle above the model ground where it left the domain, in the site
files' 20 bins of 1 km (midpoints 500, 1500, ... 19,500 m; D24). The dispersion domain top is 20 km
above ground (D21), so every side exit falls in a bin; the writer fails if one does not.

- Same reference and bins as the site files (NAME_processing reads `Z (m agl)`), and what NAME's
  fate records give directly.
- Height above sea level, or pressure, can be recreated downstream: each exit is binned by its edge
  cell, so adding that cell's terrain height gives its altitude (to within the terrain variation
  inside a cell). The file does not carry terrain heights: the inversion uses its own orography.

## Storage

Measured on our unrounded column footprints for 1-7 January 2013 (160 soundings, 2019 NAME output;
`srr` only, zlib level 4 with shuffle, one footprint per chunk):

| Encoding | kB per sounding | Column-total error | Largest cell error | Non-zero cells set to 0 |
| --- | --- | --- | --- | --- |
| float32, lossless | 70.6 | 0 | 6e-8 | 0 |
| legacy and site files: `least_significant_digit=5` (multiples of 2^-17) | 21.4 | -0.57% mean, up to -1.5% | 100% | 0.44% |
| **float32, BitRound 10 significant bits** | **41.5** | at most 4e-5 | 0.05% | 0 |
| float32, BitRound 8 bits | 35.4 | at most 1.4e-4 | 0.2% | 0 |

- **`srr`, `particle_locations_*` and `mean_age_particles_*`:** float32, BitRound to 10 significant
  bits (netCDF quantize; the file records `_QuantizeBitRoundNumberOfSignificantBits`), zlib level 4
  with shuffle, one footprint per chunk. Each value keeps a relative precision of 0.05%, far below
  the particle noise; rounding errors are unbiased, so totals are unaffected.
- **Everything else** (met, weights, release layers, coordinates, time): stored exactly; they are
  small.
- Readers need nothing special: the values are ordinary float32 with zlib. Writing needs netCDF-C
  4.9 or later (we have 4.10.1).
- Not used: zstd (about 9% smaller than zlib, but readers would need netCDF-C 4.9+ with the HDF5
  plugins; the site files were written with 4.8.1); integer packing with scale factors (fixed
  absolute precision, the same problem as `least_significant_digit`).
- Expected size: about 70 MB per month (1,285 soundings: `srr` about 53 MB, boundaries and ages
  15-20 MB), against about 115 MB lossless. The legacy month was 25 MB, but with 678 binned entries
  and the rounding bias.
- The site files carry the same rounding: in both examples the median non-zero `srr` is 6 rounding
  steps, and 10% of non-zero cells hold exactly one step, so the site footprints are probably biased
  low in the same way.

## Global attributes and provenance

The site files' attributes, plus what is needed to trace a file back to its inputs and code and to
check it matches its observation file:

| Attribute | Example |
| --- | --- |
| `species`, `model`, `met_model`, `domain`, `site` | inert, NAME, UMG, SOUTHAMERICA, GOSAT1v10-BRAZIL-column |
| `model_version` | NAME III (version 8.7.1), ACRG fate-record patch (X20) |
| `lpdm_native_output_unit` | ppm s |
| `author`, `publication_acknowledgement` | as the site files |
| `product`, `product_version`, `l2_files` | UoL GOSAT OCPR, 10.0, the L2 file names used |
| `max_level` | 16 (levels modelled by NAME, D21; the inversion checks it against the observation file) |
| `met` | UMG_Mk6PT |
| `max_travel_hours`, `integration_window` | 720, all_releases (X18) |
| `particles_per_release` | e.g. 30/s for levels 1-13, 10/s level 14, 3/s levels 15-16, 60 s (D22) |
| `random_seed` | Fixed |
| `deviations_applied` | e.g. X7, X18, X21 (IDs from docs/deviations.md) |
| `software`, `software_commit` | name-satellite, git commit |
| `created`, `Conventions` | ISO time, CF-1.8 |

## Differences from the site files

| Site file | Satellite file | Why |
| --- | --- | --- |
| dimension `time` | dimensions `sounding_id` and `footprint`, with `time` and `footprint_index` variables | times repeat (D18); footprints can later be shared |
| four faces sum to 1 | faces + top + inside sum to the column weight | column weighting, as the legacy satellite file; nothing dropped (D13) |
| `release_height` | `release_pressure`, `release_thickness`, `column_weight` per level | a column has 16 release layers |
| `least_significant_digit` rounding | BitRound, 10 significant bits | X7 |

## Differences from the legacy satellite file

| Legacy (2019) | New | Change |
| --- | --- | --- |
| dimension `time` | `sounding_id`, `footprint` | D18, D23 |
| `fp(lat, lon, time)` | `srr(footprint, latitude, longitude)` | site-file name and order; no biased rounding (X7); window covers every sounding (X18) |
| `lat`, `lon` | `latitude`, `longitude` | site-file names |
| `particle_locations_*(height, edge, time)` | `particle_locations_*(footprint, height, edge)` | every particle once (D13, X6, X14); heights above ground in 1 km bins, as now (D24) |
| none | `mean_age_particles_*`, `particle_locations_top`, `particle_fraction_inside` | D13 |
| `temperature` (°C labelled K), `pressure`, `wind_direction`, `PBLH` | `air_temperature` (K), `air_pressure`, `wind_from_direction`, `atmosphere_boundary_layer_thickness` | site-file names; X17 |
| none | `column_weight`, `release_pressure`, `release_thickness`, `release_group`, `surface_pressure` | how each column was built |
| `lev` = "column" | removed | always a column |
| attributes `max_level`, `author`, `created` | site attributes plus provenance | |

## Changes needed in the inversion code

1. Join observations and footprints exactly on `sounding_id` (both files carry it), and fail on any
   unmatched ID. Today the satellite path joins them by time, forward-filling within 60 s
   (`acrg/name/name.py` ~1013-1037, `combine_datasets` ~723), which can silently attach the wrong
   footprint now that times repeat (X22).
2. Look footprints up through `footprint_index` (one line, above).
3. Use `sounding_id` wherever `time` is the observation dimension (e.g. the `H` reshape ~1305); keep
   `time` for selecting fluxes by month.
4. Boundary conditions (~1208): the curtains go on the file's `height` bins (m above ground, as for
   the site files); decide how `particle_locations_top` and `particle_fraction_inside` enter (D13).
5. Keep the `max_level` check between observation and footprint files.

## Open questions for the group

- [ ] Monthly files as now, or one file per day (simpler to rerun a single day)?
- [ ] File name: the exact `site` part (product, version, region, "column").
- [ ] `particle_locations_top` on the full latitude × longitude grid is as large as `srr`; is a
  coarser grid, or a few latitude bands, enough for the stratospheric boundary term?
- [ ] Share footprints within release groups (through `footprint_index`) once the error is
  measured, or always store one per sounding?
- [ ] Store per-level footprints as well (about 8 times larger), e.g. to test other averaging
  kernels later?
- [ ] Time-resolved footprints: how far back (more than the site files' 24 h), in what bins, and
  whether `srr_residual` is needed (it equals `srr` minus the time-resolved sum, and is a third of
  the co2 example file). To settle after the rest of the file.
- [ ] Who updates the inversion code, and against which branch?
