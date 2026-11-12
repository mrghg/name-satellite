# How NAME 8.7.1 behaves, as observed in our runs

What we have learned from running NAME 8.7.1 on JASMIN: its inputs, outputs, error messages and
our own measurements. NAME itself is Met Office software (Crown copyright) used under licence; this
repository contains none of its code, and notes on its source are kept privately, outside the
repository.

## Particle fate output ("Fate Info? = Yes")
- One comma-separated record per particle when it stops (leaves the domain, or reaches the maximum
  travel time): `Source, Travel Time (hr), X (Lat-Long), Y (Lat-Long), Z (<units>)`, where Source is
  the integer source number. With `Met? = Yes` the met at that point follows (U/V/W Ambient, ...,
  Temperature, Pressure, ...) and each line ends with a comma. Pressure was present for every record.
- The record has no particle identifier and no reason for stopping, so these must be inferred from
  position and travel time.
- In the legacy ACRG 7.2 runs every particle appears exactly once (60 × rate + 1 per source).
- Particles released by a 60-s source are not always 60 × rate + 1: NAME 8.7.1 released 1,801 at
  30/s and 181 at 3/s, but 600 (not 601) at 10/s (smoke test 2026-11-11, job 67020916, all 22
  sources at 10/s). So the number released must be counted from the fate records (each particle
  appears once after `drop_repeated_fate_records`), not computed from the rate.
  NAME 8.7.1 writes extra records (X20):
  - 24-h benchmark: only for releases above the 19 km agl domain top (22 of 425 sources);
  - 30-day stage 1 runs: 406 of 425 sources, +12,348 records (2.2%), mostly lateral exits;
  - minimal example (2 sources, 1 thread, 30 days; `dev/scripts/prepare_fate_test.py`,
    `analyse_fate_test.py`): 3,717 records for 3,602 particles (+3.2%); 83 near-duplicate pairs
    (same source, one 10-s time step apart, positions one step apart); the excess clusters just
    before whole clock hours; records are written only at some phases of the 6-min sync cycle.
    Not a threading effect. Reported to the NAME team.
  - With the ACRG-patched build (`met.PATCHED_NAME_DIR`, outside the repository), each record also
    has `Particle Index` (1 to the number released, per source) and `Fate Code` (8 = left the
    domain, 16 = reached the maximum travel time) after Z. The minimal example then gives records
    identical to the stock build, every particle present, and 115 particles written exactly twice:
    10 s apart, both "left the domain", just before a whole hour of travel. Keeping each
    particle's first record removes them exactly (`name_output.drop_repeated_fate_records`).
  - 30-day test day with the patched build (1 Jan 2013, GOSAT-BRAZIL, 352 sources, 6-min sync
    time; procdev_20261112): 545,179 records; 8,929 particles (1.7%) written exactly twice, both
    "left the domain", but **180 s apart in every case** (not 10 s), and only 5% of pairs straddle
    a whole hour of travel. The first record is always outside the domain; the second is back
    inside it in 4.5% of pairs. So the gap is one particle time step, which varies; the first
    record is the exit. `drop_repeated_fate_records` allows gaps up to the sync time.
  - Fate codes are bit flags: one particle in that run has code 24 (= 8 + 16), at exactly the
    maximum travel time and beyond the east edge: it left on its last step. Treated as an exit.
  - In the same run 31 "left the domain" records lie exactly on the east edge as printed (positions
    have 5 decimal places), so beyond it by less than the printed precision. Counted as exits
    through that edge (`processing.classify_fates`).

## Limits met in our runs
- Footprint field requests: at most 3,000 per output time, so at most 3,000 sources per run.
- Fields output lines: a fields file with 2,720 field columns made the run stop at the END with
  "FATAL ERROR: too long an output line required for output group" (after 8 h 21 min, job
  66323622); 2,000 columns per file works. `name_input` writes groups of at most 2,000.
  Lesson: test a new run layout with a 1-hour maximum travel time on the debug queue first.
- Release-point met (`release_met.py`): at most 200 time grids, 12 arrays, 86 values per array and
  5 location sets per run; output groups default to 100 (`Max # Field Output Groups`).
- Threads: NAME will not use more threads than OpenMP offers, so the job's CPUs, `OMP_NUM_THREADS`
  and the input file's `Threads` must be equal (`lotus.py`).

## Output formats
- Output Format letters: I = indices, A = align, Z = zero lines suppressed, 2 = NAME II text.
- Fields-file headers (NAME III format) print "X/Y grid origin" equal to the input's `X Min` /
  `Y Min`, which is also the printed coordinate of cell 1: the origin is the centre of the first
  cell (stage 1 8.7.1 and legacy 7.2 files: −91.32900). NAME II-format headers were not checked
  (relevant to X3).
- Point output on a `Locations:` set (unstructured horizontal grid) and on a vertical `Z-Array`
  works. A time grid given as a `T-Array` cannot be used for field output: the run aborts with
  "UNEXPECTED FATAL ERROR in MultiplyShortTime". Use one regular time grid (nT = 1) per time.
- Point-output files (`<group>_C1.txt`): fields-file layout with 15 header rows (request, quantity,
  units, ..., "Z = <value> <units>"), preliminary columns `T Index, Location Index, T, Location`,
  one row per location, 7 significant figures.
- A pressure (Pa) point below the model ground is reported at the ground: its "Pressure (Pa)" is
  exactly the surface pressure, and temperature etc. are the surface values.
- A met-only run works with a single dummy source of one particle (1 h maximum travel time) and
  "Pressure (Pa)" requested at 0 m agl on a `Locations:` set, one time grid per sounding: 46
  soundings of 1 Jan 2013 took 15 s on 1 thread (`surface_pressure.py`). Its surface pressures are
  those the main run reports at the same points and times.
- NetCDF field output needs the NetCDF executable (`nameiii_eccodes_NetCDF_64bit_par.exe` in
  8.7.1); Mk12 met definitions exist only in 8.7.1.

## Threading benchmark (2026-11-09; 1 Jan 2013, 425 releases, 24 h, Mk6, debug queue)
| Threads | Total | Particle loop | Other | Speed-up | Efficiency |
|---|---|---|---|---|---|
| 1 | 813 s | 639 s | 174 s | 1× | 100% |
| 4 | 253 s | 162 s | 91 s | 3.2× | 80% |
| 8 | 179 s | 92 s | 87 s | 4.5× | 57% |

The particle loop scales almost linearly; the rest (met reading, output) barely changes. With
`Random Seed = Fixed`, fields and fate files are bit-identical for 1, 4 and 8 threads. NAME's timer
summary is in `<input>Timing.txt`; `sacct` TotalCPU shows 0 under mpirun.

## Runs that look successful but are not
NAME can exit 0 while writing errors to `<input>Error.txt` (e.g. particles lost when two jobs restored
met into the same folder at once). Judge a run by its error file; our job scripts fail if it is not
empty (`lotus.py`).
