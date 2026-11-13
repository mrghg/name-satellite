# Future method: scoping (branch `future-method`)

The project has three levels:

1. **Reproduce** the legacy satellite pipeline (done: stages 1–2, docs/plan.md).
2. **Improve** it and fix its bugs, keeping the footprint file as close as possible to what the
   inversion code reads now (main branch: D-decisions, X-deviations, the D23 file).
3. **A new, future-proof method**: this note. Nothing here changes level 2. Work stays on this
   branch until the group decides.

## The idea

Today a run's dispersion domain is fixed per site: particles stop when they cross its sides or
top, the face they cross gives the boundary-condition (BC) sensitivity ("curtains"), and the
footprint is written on one grid covering that domain. Domains are hard to change once runs
exist, and every site needs one chosen in advance (open question 11).

The FLEXPART-style alternative:

* **No side boundaries.** Particles are followed for the whole back-run wherever they go. The top
  stays at 20,000 m agl (D21), or could also go.
* **BC from 3D end positions.** Every particle's position (lon, lat, height) and time at the end of
  the back-run, binned on a coarse 3D grid, is the sensitivity to a global 3D CH4 field at that
  time (e.g. CAMS or our own global product), as FLEXPART's sensitivity to initial conditions.
* **Nested output grids.** Footprints on several grids whose cell sizes are multiples of each other
  and whose edges align on whole degrees (e.g. 0.1°, 0.5°, 2.5°), so finer cells sum exactly
  into coarser ones. **Both the fine and the medium grid follow each release group** (0.1°, ±5°;
  0.5°, ±25°, about 2,800 km); only the coarse grid is shared and global. Nothing is defined per
  site. (The first tests used a fixed 0.5° box over South America for simplicity.)
* **Tiles for storage and the inversion.** Processing cuts the grids into fixed geographic tiles
  and stores only those a sounding touches; the inversion chooses its region and resolution
  afterwards.

## What we already know (1 Jan 2013 test day, Nov 2026)

* **The met files are global**, and our runs declare them over their whole domain. Peak memory of
  the full day (352 sources, 4 threads) was 5.4 GB. We expected an unbounded domain to cost no
  extra memory; the smoke tests showed a fixed extra ~5.6 GB (Results).
* **Restoring met is a disk question, not a memory one.** An unbounded 30-day run needs the same
  met period as today's runs (553 GB restored for January 2013). Later periods may not need a
  restore at all.
* **Most of today's run time is a fixed cost.** A run of only 2 releases over 30 days took 31–37
  min on 1 thread; the full day (352 releases, 4 threads) took 82 min. Reading 30 days of global
  met costs about half an hour whatever the number of particles.
* **Particles stop early today.** Mean travel time at the fate record was 159 h (median 92 h, 90th
  percentile 372 h) of the 720 h maximum; 1.2% of the column was still inside at 720 h. Without
  side boundaries every particle runs 720 h: ~4.5× the particle-hours. With the fixed cost above,
  a full day would take roughly 3–4 h on 4 threads rather than 82 min (to be measured).
* **End positions are already written.** The fate records of the patched NAME (X20) give each
  particle's position when it stops, including at the maximum travel time.
* **Each output request names its own grid.** In NAME's input every field request names one
  source and one horizontal grid, so each release group can have its own small fine grid. If
  memory follows the grids each source is written on, the fine region follows each group without
  any spatial batching of releases.
* **The met resolution improves through the record** (global UM ~0.35° × 0.23° in 2013, ~0.1° in
  recent years). A footprint grid fixed at the finest level for the whole record gives the
  inversion one consistent grid; early periods are oversampled.

## Using the nested grids in an inversion

Two kinds of grid must not be confused:

* the **inversion's grid**: where the inversion solves for fluxes at 0.1°, where at 0.5° and
  where at 2.5°, chosen once for all soundings (made of fixed tiles);
* each **sounding's output grids**: 0.1° within ±5° and 0.5° within ±25° of its release (shared
  by its release group), and the global 2.5° grid.

To fill one sounding's row of H, take each flux cell of the inversion and use **the finest grid
the sounding has there; if that is coarser than the flux cell, spread its value evenly** over
the flux cells it contains:

1. a 0.1° flux cell inside the sounding's 0.1° grid: its 0.1° footprint;
2. a 0.1° flux cell outside it: the sounding's 0.5° value, spread over the 25 fine cells (the
   footprint is assumed smooth within a 0.5° cell; this only happens ≳450 km from the
   sounding, while everything within 300 km is always on its fine grid);
3. a 0.5° flux cell inside the sounding's 0.5° grid: its 0.5° footprint (equal to the sum of its
   0.1° cells, to ~0.1%);
4. a 0.5° flux cell outside it: the sounding's 2.5° value, spread over the 25 cells (≳2,500 km
   from the sounding);
5. a 2.5° flux cell: its 2.5° footprint.

The inversion's region can therefore be anywhere: each sounding brings the grids it has. The
file records, per sounding, the offsets of its 0.1° and 0.5° grids in global grids of those
resolutions.

This needs every level's cell edges on the coarser level's lines, so no cell straddles two
coarser ones. All grids share the origin (−180°, −90°) and 2.5 = 5 × 0.5 = 25 × 0.1, and **both
release-centred boxes snap outward to the 2.5° lines**: each covers whole 2.5° cells (the ±5° box
becomes 10–12.5° wide, at most 1.6× the cells). A 2.5° cell is then the natural storage tile, and
every coarse cell is either wholly inside a finer box or wholly outside it. (The first tests'
fine boxes were on 0.1° lines only.)

## Batching releases

If per-group fine grids are cheap, runs should be batched by **time**, not space: the shared cost
of a run is reading the met period, which every release in that period uses wherever it is. A
day (per region, or globally) per run, sized by run time and threads.

If they are not cheap, the fallback is spatial batching: releases grouped in predefined release
tiles, one fine grid per run covering the tile plus a margin. The margin (≈5°: particles travel
roughly 400–900 km a day near the surface) sets the fine region's size more than the spread of
the releases does.

## Consequences for the inversion (to discuss with the group)

1. **Emissions outside the region.** With no side boundaries, particles keep collecting emissions
   for the whole back-run wherever they are. Those enter through the coarse global grid and need a
   prior (fixed, or solved for coarsely). Today they are folded into the curtains.
2. **BC term.** H_bc becomes a sum over end positions of a 3D CH4 field 30 days before the
   sounding, instead of exit curtains × curtain fields. The inversion can still scale BC
   "regions" by grouping end positions (a few regions × height bands), similar to the four faces
   today, but the input format and code change.
3. **Back-run length becomes a free choice.** With BC from a 3D field, 30 days is not required;
   FLEXPART inversions often use ~10 days (~1.5× today's particle-hours, less fixed met cost).
   Shorter runs put more weight on the global model's transport and less on NAME's.
4. **Mixed resolution.** Where a sounding has only a coarser grid in a fine tile, the inversion
   spreads that value evenly over the fine cells (the far-field approximation).
5. **Particle noise at 0.1°.** Each fine cell gets 5–8× fewer particles than today's grid; the seed
   pair tests (Phase 4, D20 test) show single-cell noise grows, so fine grids may need more
   particles per release near the sounding.

## Tests

`dev/future/prepare_unbounded_test.py` prepares small runs from an existing prepared day: three
release groups of 1 Jan 2013 (48 sources, ~73,000 particles, ~14% of the day), patched NAME, 4
threads. Variants: `bounded` (today's set-up), `shared_region` and `shared_global` (one 0.1° grid
for all sources), `per_group` (±5° 0.1° grid per group), `per_group_240h`; the unbounded variants
also write the 0.5° regional and 2.5° global grids.

* **Smoke** (`--smoke`, 1 h travel, debug queue): does NAME accept the unbounded domain and many
  grids; memory per variant (set by the met and the grids from the start).
* **Full** (720 h and 240 h): run time, output size, end positions (all particles at 720 h?),
  behaviour at the poles and the dateline.

## Results (13 Nov 2026)

Runs in `/work/scratch-pw5/mattrigby/name-satellite/` (`future_smoke4_20261113`,
`future_full_20261113`); figures and `summary.json` from `analyse_unbounded_test.py` in
`future_figures_20261113`. 48 releases, 73,125 particles per run, 4 threads.

NAME syntax: an unbounded domain is written with `H Unbounded? = Yes` and the H-Coord and all
horizontal bounds left empty; NAME stops with "H-Coord, Set of Locations, and Location of Centre
should not be given for unbounded domains" otherwise.

**Memory (F1).** Peak resident memory of NAME, sampled inside the job (VmHWM; Slurm misses short
runs and timing `mpirun` measures mpirun):

| Run (1 h travel) | Peak memory | Run time |
|---|---|---|
| bounded, today's grid | 3.8 GB | 0.2 min |
| unbounded, 0.1° grid per group (±5°) + 0.5° region + 2.5° globe | 9.4 GB | 0.4 min |
| unbounded, one 0.1° grid over the region + 0.5° + 2.5° | 10.1 GB | 1.1 min |
| unbounded, one 0.1° grid over the globe + 0.5° + 2.5° | 15.6 GB | 6.9 min |

* Dropping the side boundaries costs a fixed ~5.6 GB from the start, presumably because NAME then
  works with the met over the whole globe. It does not grow as particles spread: the 720 h run
  peaked at 9.44 GB.
* Each output grid costs ~20 bytes per cell per release (0.7 GB for 720,000 cells × 48 releases;
  6.3 GB for 6.5 M cells × 48). Per-group fine grids are almost free; a shared global 0.1° grid
  is not an option (memory, and 18× slower even for 1 h).
* So releases can be batched by time (F1 answered): the fine region follows each group.

**Run time and output (F2).**

| Run (full) | Run time | Peak memory | Output |
|---|---|---|---|
| bounded, 720 h | 54 min | 5.0 GB | 33 MB |
| unbounded, 240 h | 31 min | 9.4 GB | 20 MB |
| unbounded, 720 h | 92 min | 9.4 GB | 57 MB |

1.7× today's time at 720 h, not the 4.5× of the particle-hours: the fixed cost of reading 30
days of met is large. A whole day of 352 releases, if the fixed cost is as today's (~50 min,
from 54 min for 48 releases and 82 min for 352), would take roughly 5–6 h at 720 h and ~2 h at
240 h. To be measured with a whole-day run.

**End positions.** Today 98.6% of particles leave the domain (mean 155 h); 1.4% are inside at
720 h. Unbounded: 99.9% reach 240 h (0.14% leave through the 20 km top), 99.4% reach 720 h (0.6%
through the top). At 720 h the end positions cover the globe; 35% of the column weight ends
north of the equator. Every particle has one fate record after X20 de-duplication.

**Inside and outside the region (F4, F5).** Column footprint totals, (mol/mol)/(mol m⁻² s⁻¹),
weights of each group's own sounding:

| Group | Today | 240 h inside / outside | 720 h inside / outside |
|---|---|---|---|
| NE coast of Brazil | 0.48 | 0.39 / 0.52 (57% outside) | 1.03 / 4.01 (80%) |
| Amazon | 2.23 | 2.43 / 1.08 (31%) | 3.02 / 5.60 (65%) |
| northern Argentina | 0.86 | 1.02 / 1.24 (55%) | 1.67 / 5.41 (76%) |

* Over 30 days most of the sensitivity is to fluxes outside South America; today that is hidden
  in the curtains. Without boundaries the outside fluxes need a prior, or a shorter back-run.
* Inside the region the 720 h footprint totals 1.4–2.2× today's: particles that leave and come
  back pick up regional emissions that today's runs give to the boundary term. **All of the
  extra is far from the sounding** (below).
* Caveat: the test region (medium grid) extends ~0.5° beyond today's domain, and its edges are not
  on 2.5° lines, so coarse cells along the edge are assigned by their centre.

**By distance from the release** (column footprint inside the region; `by_distance` in
`summary.json`, `fu6_by_distance.png`). Bounded runs on today's grid, unbounded on the 0.5° grid,
so the smallest band is 300 km (smaller bands would differ by binning alone):

| Group | < 300 km: today / 240 h / 720 h | 300–1,000 km | 1,000–3,000 km | > 3,000 km |
|---|---|---|---|---|
| NE coast of Brazil | 0.089 / 0.085 / 0.091 | 0.18 / 0.15 / 0.20 | 0.10 / 0.13 / 0.42 | 0.10 / 0.03 / 0.32 |
| Amazon | 0.164 / 0.166 / 0.168 | 0.35 / 0.33 / 0.35 | 1.19 / 1.07 / 1.29 | 0.53 / 0.86 / 1.21 |
| northern Argentina | 0.164 / 0.162 / 0.168 | 0.14 / 0.15 / 0.18 | 0.55 / 0.64 / 0.99 | 0.01 / 0.08 / 0.33 |

Within 300 km, 240 h is within 4% of today and 720 h within 2–7% of 240 h. The extra of 720 h is
beyond 1,000 km: sensitivity to fluxes far away and weeks earlier, which a shorter run gives to
the boundary term.

**Proposed default: 10-day (240 h) back-runs** (Matt, 13 Nov 2026), pending F9: same near field
as 30 days, cheaper than today's runs (31 vs 54 min here), and 31–57% rather than 65–80% of the
footprint outside the region. The cost is that the boundary term comes from a global CH4 field
10 days before the sounding, so its accuracy matters more.

**Nesting.** Summing 5 × 5 fine cells (area-weighted) reproduces the 0.5° cell to 0.06% (weighted
over cells; median 0.09%, max 1.3%, per release, Amazon group): close, not bit-exact, probably
because the conversion to mixing ratio is per cell. The test's fine boxes were on 0.1° lines only;
a real scheme aligns every level.

## Open questions

* F1. ~~Does NAME's memory follow the grids each source is written on?~~ Yes (Results).
* F2. Run time of a whole unbounded day; how much is fixed met cost (estimate 5–6 h at 720 h).
* F3. Which global 3D CH4 field for the BC, and how the inversion scales it.
* F4. Prior and treatment of emissions outside the inversion region.
* F5. Back-run length: 240 h proposed (Results); confirm with F9.
* F6. Particles per release needed for 0.1° footprints near the sounding.
* F7. Tile scheme (sizes per level) and the file format the inversion would read.
* F10. Size of the release-centred 0.5° grid (±25° proposed): measure memory, run time and how
  much of a 10-day footprint falls outside it; compare with one shared global 0.5° grid.
* F8. Why nested grids differ by ~0.1% (mixing-ratio conversion per cell?), and whether to
  output a quantity that adds exactly.
* F9. Boundary term: sample a global CH4 field (e.g. CAMS) at the end positions and compare
  with today's curtain term for the same soundings.
