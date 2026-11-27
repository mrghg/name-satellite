# Sizing future-method runs (draft, 27 Nov 2026)

How many soundings to put in one NAME run, for the fewest core-hours. Built from the January and
August 2013 test runs (10-day back-runs, 4 threads, `future-method` set-up); good to about ±25%,
to be refined as more runs are made. The NAME behaviour behind it is in
`dev/notes/name-behaviour.md` ("Run time and field requests").

## What a run costs

NAME's run time, 4 threads, 10-day back-runs:

    T ≈ F + a·P + k·P·N

| Term | Value | What it is |
|---|---|---|
| F | 3,000–4,500 s per run (take 3,500) | Reading the met: ~134 GB for 11 days of global UM met, Mk6 and Mk7 alike. Varies with read speed (4–14 GB/min); slower when several runs read at once. |
| a·P | a ≈ 4,400 s per million particles | Moving the particles. P = particles released in the run. |
| k·P·N | k ≈ 7 s per field request per million particles | Measured: run time grows with particles × field requests in the run, whichever source or species each request is for. N = field requests in the run. |

The terms come from: the 30 January days (one run each), 7 January split into four runs, the
same 240 releases with 720 and with 126 field requests (particle loop 3,444 s and 1,942 s), and
two August (Mk7) runs. For one run per January day the formula gives 83 h of NAME time; the runs
took 74 h.

## Layouts

* **species** (the default from 27 Nov 2026, D28): one species per sounding, three field
  requests per sounding (0.1° box, 0.5° box, global 2.5° grid). N = 3 × soundings.
* **sources** (used for the January 2013 month): three requests per release, 16 releases per
  release group. N = 3 × releases, so 6–11× larger; the k·P·N term then dominates big runs
  (about 70% of a 944-release run).

Both give the same column footprints (7 Jan 2013 part 4, 42 soundings: within 3 × 10⁻⁶ of each
field's maximum, end positions identical).

## Best run size (species layout)

With S soundings in a run and p particles per sounding (January 2013: p ≈ 12,800), the cost per
sounding is F/S + a·p + 3k·p·S, smallest at

    S* = √(F / (3·k·p)) ≈ 115 soundings

and flat around it: half or twice S* costs about 13% more per sounding. So:

* **Aim for about 100–200 soundings per run**, taken in time order (neighbouring soundings share
  met and release groups).
* **Combine consecutive days that are small** (January averaged 42 soundings a day): a run
  spanning d days reads met for 10 + d days, about 10% more per extra day, which is small next to
  F saved by running fewer runs. For January 2013: one run per day ≈ 53 h of NAME time, two-day
  runs ≈ 44 h, three-day runs ≈ 44 h (estimates).
* **Split days with more than ~200 soundings.**
* Limits per run: at most 216 soundings (species particles can carry; not yet reached), 3,000
  field requests, 7 M particles (`MAX_PARTICLES`), 24 h on the `high` queue. Memory was
  ~10.5 GB in every test run (42 species at most so far).

## Still being measured (Nov 2026)

1 and 2 threads (core-hours per run; 8 threads were 1.45× faster than 4, so 73% per core), four
identical runs started together (read contention), 3–5 Jan as two runs (met for 2–3 days per
run), and the global grid only (one request per release). Update F, a, k and S* from these and
from the first production months.
