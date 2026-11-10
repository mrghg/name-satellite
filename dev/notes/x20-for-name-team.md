# Extra particle-fate records in NAME 8.7.1 (X20): minimal example for the NAME team

Written 2026-11-10 for Matt Rigby to send to the Met Office NAME team.

## Summary
In backward runs with stock NAME 8.7.1, the particle "Fate Info" output contains **more records than
particles released**: about 3% more over 30 days. In the legacy ACRG NAME 7.2 runs every particle
appeared exactly once (60 × rate + 1 records per source). The extra records look like the same
particle written twice, one model time step apart. They are concentrated around whole clock hours.
The fields are not affected; we use the fate records to compute boundary sensitivities (where each
particle left the domain), so we need each particle once.

## Minimal example
- NAME 8.7.1 (stock), **one thread**, LOTUS2 (JASMIN), job 66602864, run time 31 min.
- Backward run, Global Mk6 met; 2 sources (one GOSAT sounding, two 60 s releases at 30 particles/s:
  1,801 particles each, 3,602 in total), 30-day maximum travel time, sync time 6 min.
- Output: `Sets of Particle/Puff Information` with `Fate Info? = Yes`, `Met? = Yes`, Z in m agl.
- Files: `/work/scratch-pw5/mattrigby/name-satellite/fate_test_20261110/` (`input.txt`,
  `output/Particles_C1.txt`, `release_plan.csv`); made by `dev/scripts/prepare_fate_test.py`,
  analysed by `dev/scripts/analyse_fate_test.py`.

## What we see
- **3,717 records for 3,602 particles** (+115, 3.2%): 1,860 and 1,857 for the two sources.
- Fates: 3,660 lateral exits (within ~0.002° of a domain edge), 57 at the maximum travel time, none
  through the top, none unexplained.
- **Near-duplicates:** 83 pairs of records from the same source, one 10-s time step apart in travel
  time (Δ = 0.0028 h), positions about one step's motion apart (< 0.01°, < 5 m). No exact
  duplicate rows. Example (Source, travel time h, lon, lat, z m agl):
  `1, 44.99530, -19.44744, -14.91287, 709.34` and `1, 44.99807, -19.45003, -14.91460, 705.20`.
- **Timing:** by clock time (release time − travel time), records appear only at some phases of the
  6-minute sync cycle, and the counts rise just before each whole hour. Records per 2-minute bin,
  00-02 … 58-60 min past the hour:
  `175, 0, 228, 126, 0, 226, 107, 0, 235, 120, 0, 235, 113, 0, 228, 122, 0, 224, 111, 0, 230, 103,
  0, 253, 123, 0, 277, 126, 0, 355`.
  So the near-duplicate pairs sit just before whole hours (e.g. travel times 44.995 h, 50.997 h,
  53.99 h, 77.99 h).
- Not a threading effect (this run is single-threaded); in our 4-thread 30-day runs (425 sources)
  the excess was 2.2%. In 24-h runs the excess appeared only for releases above the domain top.

## Our reading of the code (please correct)
Fate records are written in `CalcParticleResults` (Case.F90 ~6098-6119) for marked particles
(`Particle%Marked /= 0`) whenever the particle's time centre crosses a requirement time
(`ReqInfoT`, and the travel-time list `ReqInfoS`, Case.F90 ~4150-4209). Marked particles are only
inactivated at the end of the sync step (Particle.F90 ~127-140). We suspect a marked particle can
take one more step and cross a second requirement time (e.g. at a whole hour) before it is
removed, and is written again.

## Questions
1. Is this expected? Can each particle be written exactly once?
2. Could the fate output include a per-particle index (as the non-fate particle output has
   "PP Index"), so duplicates can be removed reliably?
3. Is there a way to get the fate code (left the domain laterally / through the top / maximum
   travel time) in the output, rather than inferring it from position and travel time?
