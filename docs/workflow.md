# Workflow

This page describes the end-to-end workflow in plain terms. It describes the **target** design;
steps not yet built are marked *(planned)*. Each step names the module that does it, so you can go
from here to the code. Each module's docstring describes the same step in more detail.

```
 Satellite L2 files ─► 1. Read soundings ─► 2. Plan releases ─► 3. Write NAME inputs
                                                                        │
         footprint netCDF ◄─ 6. Combine into column ◄─ 5. Read NAME output ◄─ 4. Run NAME (LOTUS2)
```

## Preparing runs from current L2 data (`name_satellite/prepare.py`)
`python -m name_satellite.prepare PRODUCT SITE YYYY-MM WORK [--test]` chains steps 1-4 for a
month. For each day it downloads the L2 file (`l2_files.py`, CEDA over HTTPS; days without a file
are skipped), selects the good soundings in the site's box (`domains.py`), writes the observation
file (keyed by sounding ID) and `day.json`, and writes two NAME jobs:

1. a **met-only surface-pressure run** (`surface_pressure.py`; 1 thread, minutes): the model's
   surface pressure at each sounding's position and time, needed to scale the release layers (D19),
   and the pressure at the domain top (20 km agl), to check the modelled layers lie inside (D21);
2. the **main job**, which first runs `python -m name_satellite.prepare main-input DAYFOLDER`
   (`write_main_input`: reads the surface pressure, checks model / retrieval is within 0.9-1.1,
   scales each sounding's levels (D19), checks that the modelled levels (16 for the GOSAT products)
   lie inside the domain and carry at least 0.85 of the column weight (D21), groups the soundings
   that share releases (D20), and writes
   `surface_pressure.csv`, `release_groups.csv`, `release_plan.csv` and `input.txt`) and then runs
   NAME.

It also writes one met-restore job for the month and `submit.sh`, which chains restore → each day's
surface-pressure job → that day's main job (`afterok`: a job starts only if the one before
succeeded). Job scripts exit non-zero if NAME fails or writes to its error file. `--test` limits
travel to 1 hour on the debug queue.

## 1. Read soundings *(first version: `name_satellite/observations.py`)*
**In:** one day's University of Leicester OCPR proxy XCH4 file: GOSAT-1 v10.0 or GOSAT-2 v2.0 (ESA
CCI v7.x for legacy validation). **Out:** a table with one entry per sounding (dims `sounding`,
`level`): exposure ID, time, latitude/longitude, XCH4 and its uncertainty, quality and land/glint
flags, and per-level pressure (hPa, surface first), pressure weight, averaging kernel and a priori.

1. `read_l2` reads the file, identifies the product and version, checks the layout and that
   exposure IDs are unique, and records the time resolution (v10 and GOSAT-2 store whole minutes,
   rounded down: deviation X22).
2. Filters give masks, applied with `select`: `quality_mask` (flag 0), `bad_pressure_mask` (no
   undefined level, NaN or -9999.99: X24), `level_order_mask` (pressure decreasing from the
   surface), `mode_mask` (land or glint), `latlon_mask` (lower bound included, upper excluded).
3. Optional `legacy_bin` averages the soundings in each lat/lon bin as legacy did (quirks: X23) and
   adds `dlat`/`dlon` (footprint + spread) for the release boxes.
4. `to_obs` / `write_obs_file` give our observation file, indexed by `sounding_id` (D18), with the
   same variable names as the legacy file plus `time`, and `dlat`/`dlon`/`n_soundings` for bins.
   `releases.build_release_plan` takes it and carries `sounding_id` into every release.
   `to_legacy_obs` / `write_legacy_obs_file` give the legacy time-indexed layout, for comparisons
   with legacy files only; they refuse repeated times (X22).

Checked against the legacy `acrg.satellite.gosat` on ESA CCI v7.0 files for 1 and 7 Jan 2013
(`dev/scripts/run_legacy_gosat.py`, `compare_observations_with_legacy.py`): same soundings, values
identical, unbinned and binned. The legacy code cannot read v10 or GOSAT-2 files (X22).

## 2. Plan releases *(in progress: `name_satellite/releases.py`)*
**In:** the sounding table and a pressure ceiling. **Out:** a release plan: for each sounding and
each retrieval level below the ceiling, a vertical layer (pressure top and bottom, Pa), a
horizontal box (the footprint size), a 60 s release time and a particle rate.
Layers run between the midpoints of neighbouring retrieval levels. The lowest layer starts at the
surface. For new runs every level is first scaled by model / retrieval surface pressure
(`scale_to_model_surface`, D19), so the layers follow the model's terrain and keep their shares of
the column; the model surface pressure comes from the met-only run (`surface_pressure.py`).
Soundings in the same NAME output-grid cell and within 60 s of each other (GOSAT's repeat
exposures) then share one set of releases (`release_groups.py`, D20): mean position, time and
scaled levels; box = footprint + spread; ID = earliest member's sounding ID. Nothing is averaged
in the observations: each sounding gets its own column footprint later, from its group's
per-level footprints and its own weights. `Site.release_cell` sets the cells ("grid", a
(dlat, dlon) size, or None). Legacy
validation runs instead replace level 0 by the legacy NAME surface pressure (D16, X9). Legacy
details: `dev/notes/gosat-l2-and-preprocessing.md`.

Particle rates: new runs take each level's rate from `releases.PARTICLE_RATES` (D22: 30/s for
levels 1–13, 10/s for level 14, 3/s for 15–16 for the GOSAT products); validation runs keep the
legacy rule (30/s up to level 12, 3/s above).

Done: `release_layers_pa` (layer centre and thickness, Pa), `particles_per_second` (legacy
rate rule), `footprint_box_degrees` (10.5 km box, legacy haversine) and `build_release_plan` (one
row per release: name, sounding, level, times, box, layer, rate). They reproduce the vertical layers, positions, times and rates of all 11,526 legacy
sources for January 2013 (`tests/test_releases.py`). Still to do: pressure ceiling.

## 3. Write NAME inputs *(first version: `name_satellite/name_input.py`, `met.py`)*
**In:** the release plan, domain settings and the met version available for the dates. **Out:** a
NAME input file per day: one source and one output request per release, plus particle-fate output.

`met.met_version_for_period` picks the UM version (NAME team's 8.7.1 dates);
`required_met_files`/`missing_met_files` check the archive (all 3,542 Mk6 files for 1 Jan 2013
are present). `name_input.write_name_input(path, plan, RunSettings(...))` writes the file in
NAME 8.7.1 syntax, using the NAME team's met declarations and restore script. Checked by a short
NAME 8.7.1 run on LOTUS2.

Batching (`batching.py`): several days can share one NAME run, so the met they have in common is
read once. `combine_release_plans` joins the daily plans, numbering the sources over the whole run
(each keeps its per-day name); the integration window must then cover every release
(`integration="all_releases"`). NAME 8.7.1 allows at most 3,000 sources per run, i.e. about a
week of GOSAT-BRAZIL. After the run,
`split_batched_output` cuts the fields and particles back into days that look exactly like
one-day runs, and each day is processed as usual (step 6). Checked in stage 2 (1-hour test of all
2,720 releases); batching saves little compute, so daily runs are the default. Field requests are
written in output groups of at most 2,000 (`name_input.fields_output_group`), because NAME refuses
over-long output lines (2,720 in one group failed); `name_output.read_fields_folder` reads all the
files of a run.

## 4. Run NAME *(job scripts: `name_satellite/lotus.py`)*
Each day is a SLURM job on LOTUS2. `lotus.write_job_script` writes the job with the NAME team's
environment; one thread count sets the job's CPUs, `OMP_NUM_THREADS` and the input's `Threads`. NAME runs backwards in time for up to 30 days. It writes, for
each release, a footprint field (time-integrated or hourly) and a record of where and when each
particle left the domain.

## 5. Read NAME output *(in progress: `name_satellite/name_output.py`, `legacy.py`)*
**In:** NAME fields and particle files. **Out:** for each release, a footprint (sensitivity of the
mole fraction at the release to surface fluxes, (mol/mol)/(mol/m²/s)), counts of particles leaving
through each boundary (for boundary conditions) and the met at the release point.
Everything is matched by release ID, never by order.

Done: `read_fields` (NAME text fields → array over source × lat × lon, keyed by source name;
missing all-zero rows restored), `read_acrg_particle_file` (legacy particle fates), and legacy
release-point met (`legacy.read_legacy_met_folder`). Tested on every day of January 2013. Stock
NAME 8.7.1 fate files: `read_name_fate_file`, turned into the legacy particle table by
`legacy_pipeline.legacy_particles_from_fate` (heights in m agl). Still to do: hourly fields.

Release-point met (`release_met.py`): the main run itself outputs, for each sounding, the met at
its position and time (level-1 temperature, pressure and wind, boundary-layer depth, and surface
pressure), instead of a separate legacy-style met run with one set of grids per release, which
NAME's grid limits don't allow for busy days. `read_release_met` reads it, and
`release_met_per_source` turns it into the per-release table the processing uses, with every
release's pressure = min(requested, model surface pressure) (X9). Limit: 199 soundings per run.

Particle fates (decision D13): every released particle ends in exactly one category. It either
left through the north, east, south or west face, left through the top, or was still inside at the
end of the run. Edge and top exits are stored as histograms of *fractions of the particles released*;
for particles still inside only their fraction is stored. All categories together sum to exactly 1. The inversion multiplies these by a global model's mole
fractions to get the background, and decides how to treat top exits and particles still inside.

## 6. Footprints, boundary histograms and column combination *(legacy method done)*
**In:** the per-level footprints and the sounding table. **Out:** one footprint per sounding,
Σ over levels of (pressure weight × averaging kernel × level footprint), written in the netCDF
format the inversion reads (legacy `fp` format and/or `srr`), with `max_level` recorded.

Modules, in order:
- `footprints.py` (6a): `grid_cell_area_m2`, `footprint_from_mixing_ratio`: field [ppm s] × 1e-6 ×
  cell area / release duration → (mol/mol)/(mol/m²/s).
- `boundary_histograms.py` (6b): `legacy_boundary_histograms` (variants "2019" and "2025"; the
  D13 method is still to do).
- `column.py` (6c): `nearest_release_levels` (legacy nearest-pressure rule) and `combine_levels`
  (Σ footprint × pressure weight × averaging kernel over levels below `max_level`).
- `legacy_pipeline.py`: `process_day_legacy` chains steps 5–6 exactly as legacy (reading a legacy run folder; `process_legacy_inputs` is the same processing for data already read, e.g. NAME 8.7.1 runs);
  `write_legacy_fp_file` writes the 2019 file layout and encoding.

Result: 1 Jan 2013 production footprints rebuilt **bit for bit** (every variable;
`tests/test_processing.py`). The "2025" histogram variant also matches today's legacy code exactly.

Writing: `write_legacy_fp_file` writes the legacy layout. `rounding="legacy"` reproduces the 2019
files' lossy rounding (`least_significant_digit`), which biases column totals low by 0-1.4%;
`rounding="none"` stores full float32 values and is the setting for new footprints (X7).

## Validation
Each step is first checked against the legacy pipeline's output for GOSAT-BRAZIL, January 2013
(see [reference-data.md](reference-data.md)), before any deliberate improvements
([deviations.md](deviations.md)).

NAME is a Monte Carlo model, so new runs are compared with legacy statistically:
`validation.py` gives the measures (total ratios, footprint shape differences at several grid
coarsenings, edge-histogram differences), and each is set beside the same measure for two runs
that differ only in their random seed (pure particle noise). Scripts: `dev/scripts/process_stage1.py`
(1 Jan 2013) and `process_stage2.py` (1–7 Jan, including the batching test).
