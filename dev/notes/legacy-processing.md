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

## X14: boundary-histogram edge assignment (investigated 2026-11-07)
Code: `particle_locations()`, legacy process.py ~1181 (today), acrg `acrg_name/process.py` (2019/2021).
- Edge sets: N = Lat > lat_max − dlat/2, E = Long > lon_max − dlon/2, S and W likewise, using the
  **output grid**. Each set is binned with `np.histogram2d` over (along-edge cell edges, height edges
  0..20 km). Then each (time, level) histogram is normalised by the sum over the four edges.
- 2019 (acrg d8a5e33d) → 2021 (01c44ad3): only additions (mean ages, 0–360° special case);
  N/E/S/W counts unchanged.
- **2025-05-15, NAME_processing 094d8e6** ("updated to have max of xvalues and yvalues"): the upper
  outer bin edge is stretched to `max(values)`, so particles beyond the far end are now counted in the
  last cell. The lower end is unchanged, so they are still dropped there. The benchmark was regenerated
  the next day (974c176), which is why today's output "matches the benchmark".
- 1 Jan 2013 MR run (dispersion domain = grid + 5°), out of 562,925 particles:

  | Edge | Particles | Below first cell (always dropped) | Beyond last cell (dropped 2019, kept 2025+) |
  |---|---|---|---|
  | N | 5,748 | 1,160 | 830 |
  | S | 74,461 | 5,798 | 4 |
  | E | 286,332 | 4 | 761 |
  | W | 188,306 | 5,683 | 1,170 |

  Also 8,011 corner particles are on two edges, and 16,089 are on none (still inside at max age, or
  left through the top; see X6).
- Result: same per-level totals (renormalised), but shifted edge/height fractions. Corner regions and
  edge ends are the most affected.

## Particle height reference (checked 2026-11-07)
Legacy particle-file `Ht` is **metres above ground level**:
- The ACRG-modified v7.2 build writes the particle's internal NAME height, unconverted (checked in its
  source; details in the private notes).
- Empirically, particles leaving through the top (inside the lateral domain, age 1–719 h, 1 Jan 2013 MR
  run) have median Ht 19,011 m over the Altiplano (terrain ~4 km), 19,012 m over the Atlantic and
  19,013 m over the Amazon. The domain top is 19,000 m agl.
- The processing bins heights 0–20 km in 1 km steps with no stated reference (open question 8).

## Zero surface footprints at mid levels (1 Jan 2013 MR run; to investigate with new runs)
Number of soundings (of 25) whose release at each level has an all-zero surface footprint:

| Level | 1–3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Zero | 0 | 1 | 8 | 11 | 11 | 10 | 4 | 6 | 6 | 7 | 9 | 11 | 17 | 23 | 25 |

For example, sounding 1 has nonzero totals at levels 1–3 and 9–14 but zero at 4–8. Rapid lateral
exit can explain zeros, but the non-monotonic pattern with height is worth checking: it may relate
to the "poor performance high up" anecdote, or to the 30/s → 3/s particle split (levels ≥ 13).

## Release-point met = release pressure (checked 2026-11-07)
The legacy met run samples met at each release point, with the height given in Pa, so the met
"Pressure" equals the requested release-layer centre (rtol 1e-6). The only exceptions are releases
below the model ground: 22 of 11,526 in January 2013 (levels 1–2), where NAME reports a lower
pressure (the surface). The processing's "nearest NAME level by air_pressure" therefore matches
retrieval levels to release-layer centres (except for those cases), not to an independent met
pressure. Particle counts are exactly 60 × rate + 1 per source (NAME releases at both ends of the
60 s window).
