# Legacy processing (NAME_processing)

Investigated 2026-11-06 in `~/code/NAME_processing` (commit 309e32e).
`name_processing/process.py` is about 4870 lines. Key functions (approximate line numbers):

| Function | Line | Role |
|---|---|---|
| `read_file` | 399 | NAME text fields |
| `define_grid` | 598 | grid from header |
| `read_met` | 736 | point-met files |
| `met_satellite_split` | 1136 | per-point met |
| `particle_locations` | 1181 | boundary histograms |
| `footprint_array` | 1648 | build fp |
| `footprint_concatenate` | 2281 | |
| `write_netcdf` | 2547 | |
| `satellite_vertical_profile` | 2871 | column combination |
| `define_subfolder` | 3041 | |
| `process` | 3348 | entry point |
| `find_satellite_obs_file` / `open_obs` | 4800 / 4822 | obs file `Observations/*{datestr}_*.nc` (exactly one) |

## Behaviour
- Mode: `species=None` → `MixR_files` (time-integrated); species given → `MixR_hourly`; HiTRes only for CO2
  with `user_max_hour_back > 0`.
- Units: ppm s → value × area × 1e-6 / (3600 × timeStep); g s/m³ uses molm3 = 345/R by default.
- Particle Id → (time, level) via a running counter (assumes contiguous IDs).
- Histograms (`particle_locations_n/e/s/w`, 20 heights) are normalised by edge particles only; top
  exits are excluded.
- `satellite_vertical_profile`:
  - the lat/lon check within 1° only logs;
  - the per-sounding sum of particle histograms must equal the number of levels;
  - `fp.time` is replaced by the obs time;
  - for each retrieval level 0..max_level-1, it takes the nearest NAME level by `air_pressure`
    (release-point met pressure) and weights it by `pressure_weights × xch4_averaging_kernel`;
  - checks Σ weighted == Σ(PW·AK); `lev="column"`.
- `max_level` drops retrieval levels ≥ max_level (no copying; the inversion fills them from the a priori).
- Output: `{process_dir}/{DOMAIN}/{SITE}-{height}/{species}/{SITE}-{height}_NAME_{met}_{DOMAIN}_{species}_{YYYYMM}.nc`.
  The production satellite files are named `GOSAT-BRAZIL-column_SOUTHAMERICA_201301.nc` with variable `fp`.
- `run_nc_hourly.py` (`create_layered_sparse_netcdf` ~673; grid origin +0.5·res ~340) and
  `run_nc_hourly_satellite.py` (folder `REGION_SITE[_YYYYMMDD]_SUFFIX` → `Fields_{site}_MixR_hourly_{date}.nc.bz2`).

## Suspected bugs → docs/deviations.md X3–X8, X11.

## Tests
23 satellite tests pass under `module load jaspy/3.12/v20250704`. 1 site test fails
(`test_run_nc_hourly`: 27.904488 vs 27.904209).

## Upper levels (2013-01-01, MR 2019 run, Z max 19 km)
- Levels 1–12: 1801 particles each; no age-0 deaths; 0.1–6.5% still in the domain at 720 h.
- Levels 13–17: 181 particles each.
- Level 16: 36% exit through the top.
- Level 17 (layer about 65–30 hPa, centre 47.5 hPa ≈ 20.5 km): 89% die at age 0; 99% leave through the top.
- Obs levels at 10, 1 and 0.1 hPa are above the met top; their PW·AK is about 2.4% of the total.
