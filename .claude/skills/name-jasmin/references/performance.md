# NAME 8.7.1 run time and memory (measured on LOTUS2, Nov 2026)

All runs: UM global Mk6 met, GOSAT-BRAZIL releases (16 levels per release group, 60-s releases,
3–30 particles/s per level), patched 8.7.1 unless noted.

## Threads (1 Jan 2013, 425 releases, 24 h travel, debug queue)
| Threads | Total | Particle loop | Other | Speed-up |
|---|---|---|---|---|
| 1 | 813 s | 639 s | 174 s | 1× |
| 4 | 253 s | 162 s | 91 s | 3.2× |
| 8 | 179 s | 92 s | 87 s | 4.5× |
The particle loop scales almost linearly; met reading and output barely change. With a fixed
seed, output was bit-identical for 1, 4 and 8 threads (the docs only promise near-identical). 4 threads is a good default.

## Run time: a large fixed cost per run
- 30-day back-run, 352 releases, 4 threads, bounded South America domain: **82 min**.
- Same domain, 48 releases: 54 min. Two releases, 1 thread: 31–37 min.
- So roughly 50 min of every 30-day run is fixed cost (reading 30 days of global met), whatever
  the number of particles. Batching more releases (or days) into one run amortises it.
- Without side boundaries (every particle runs the full time), 48 releases: 92 min at 720 h,
  31 min at 240 h — 1.7× and 0.57× the bounded run, not the 4.5× the particle-hours suggest.
- Particles in a bounded 30-day South America run stop after 159 h on average (median 92 h).

## What NAME's timers miss
With the default on-demand met method ("new threaded method"), the per-thread met-read timers in
`<input>Timing.txt` record no calls, so met reading is not timed. In a 240 h unbounded run of 48
sources the summary accounted for only 12% of the run (particle loop 191 s of 1,834 s); in 1-hour
tests, about 75%. Measure from outside instead: sample the NAME process's `/proc/<pid>/task/*/stat`
(CPU per thread), `/proc/<pid>/io` (`rchar`: bytes read, including from network file systems;
`read_bytes` stays ~0 on Lustre) and `VmHWM`, and count "read from file" lines in `<input>Log.txt`.
`src/name_satellite/future/sampler.py` (branch `future-method`) does this every 5 s. A 1-hour
unbounded test read ~6.5 GB of Mk6 met (`rchar`).

## A month of 10-day unbounded runs (Jan 2013, 22 days, 4 threads, Mk6)
Settings: no side boundaries, 20 km top, deep convection `New`, near-source scheme for 1 h,
0.1° and 0.5° boxes per release group plus a global 2.5° grid (`future-method` branch).
- **Run time follows the particle count**: 43 min (0.15–0.2 M particles) to 5.2 h (1.05 M).
- **Every run read ~134 GB of met** (`rchar`; ~1,170 Mk6 files for 11 days of global met), at
  6–15 GB/min while reading. Peak memory 9.5–10.3 GB.
- **Waiting on met reads is small**: periods with all threads idle while `rchar` rises are 0–5%
  of busy runs, up to ~22% on light ones (too little particle work to hide the reads). With 4
  threads, 3.0 were busy on average.
- NAME's timers: particle loop 30–65% of the total, output < 2%; the rest (25–70%) is untimed
  even with `Summary Only? = No` (met timers still record zero). That untimed time is
  multi-threaded and grows with the particle count, so it is not just reading.
- Variants, each against the same day's normal run:

| Change | Wall time |
|---|---|
| 8 threads instead of 4 | 191 vs 277 min (1.45× for 2× the cores) |
| `Prefetch? = No` | 120.1 vs 120.1 min (no effect) |
| `Deep Convection? = No` | 176 vs 207 min (convection costs ~15%) |
| Detailed timers (`Summary Only? = No`) | 119 vs 123 min (no cost) |

## Memory
- Bounded domain: 3.8 GB (1 h travel) to 5.0–5.4 GB (30 days).
- **No side boundaries adds a fixed ~5.6 GB** (9.4 GB), present from the start and not growing as
  particles spread over the globe (9.44 GB after 720 h).
- **Each output grid costs ~20 bytes per cell per field request**: one 0.1° grid over South
  America for 48 releases +0.7 GB; one global 0.1° grid +6.3 GB (and 18× slower even for 1 h).
  A small grid per release group (e.g. 0.1° ±5°, 10,000 cells) costs almost nothing.
- Particles themselves are cheap (~545,000 particles per 352 releases).

## Output size
Fields files list only non-zero cells; see outputs.md for measured sizes. A day of the bounded
30-day South America run: ~600 MB of text (fields plus particle-fate records).

## Footprint noise (for planning particle numbers)
Two runs differing only in seed: column-footprint totals differ ~3–4% per sounding; about a third
of a single sounding's footprint is "in different cells" on a 0.35° grid, falling to ~5% on 5.6°
cells. More particles in the upper levels (30/s up to level 13, 10/s at 14) cut noise variance
~37% for ~10% more compute.
