# Sizing future-method runs (draft, 28 Nov 2026)

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

## What the profiling runs showed (27 Nov 2026)

All on 7 Jan 2013 part 1 (240 releases, 0.37 M particles, 720 requests, sources layout) unless
stated. "Rest" is NAME's total minus the particle loop: mostly reading met.

| Run | Threads | NAME total (s) | Particle loop (s) | Rest (s) | Allocated core-s |
|---|---|---|---|---|---|
| alone, twice | 4 | 5,564; 6,794 | 2,571; 2,683 | 2,993; 4,111 | 22,000–27,000 |
| four started together | 4 | 5,012 ×3; 6,758 | 2,208–2,633 | 2,729–4,125 | 20,000–27,000 |
| 2 threads | 2 | 7,488 | 4,790 | 2,698 | 14,980 |
| 1 thread | 1 | 12,055 | 9,116 | 2,939 | 12,060 |
| global grid only (240 requests) | 4 | 3,098 | 1,551 | 1,547 | 12,400 |
| 3–5 Jan as two runs (256, 272 releases; met for 13 days) | 4 | 7,714; 8,355 | 3,446; 3,466 | 4,268; 4,889 | 31,000–33,000 |

* **Threads:** the particle loop scales almost linearly (4 threads 3.7× faster than 1; 8 threads
  only 1.45× faster than 4, from an earlier test), but reading the met does not speed up at all.
  So **one thread costs about 40–50% fewer core-hours** than four for the same run. The 1-thread
  job kept its core busy 87% of the time; the 4-thread jobs 78%.
* **Contention:** four runs reading met at once were no slower than one alone. The spread of the
  "rest" (2,700–4,100 s for identical runs) is the storage or the node, not each other: the one
  slow run of the four was on the node where the 2-thread and 3–5 Jan runs also ran.
* **Multi-day runs:** a run spanning 3 days spent 4,300–4,900 s outside the loop against
  2,700–4,100 s for one day: roughly the expected ~10% more per extra day of met.
* **Global grid only:** cutting 720 requests to 240 saved 1,020 s of loop: k ≈ 5 s per request per
  million particles here, and a ≈ 3,000 s per million particles. Both are within the run-to-run
  spread of the earlier fit (the four parts of 7 Jan, with the same particles and requests, took
  2,150–3,700 s in the loop). Its "rest" was unusually short (1,547 s), with the same met: read
  speed varies by more than the requests' share of it.
* **Memory:** 10.3–10.8 GB in every run (any thread count, up to 42 species).

## Best run size (species layout)

With S soundings in a run and p particles per sounding (January 2013: p ≈ 12,800), the cost per
sounding in core-seconds, with c threads, is c·F/S + a₁·p + 3·k₁·p·S, where a₁ ≈ 14,000 s and
k₁ ≈ 22 s per million particles are the single-thread values (about 3.7× the 4-thread a and k).
It is smallest at

    S* = √(c·F / (3·k₁·p))   ≈ 60 soundings with 1 thread, ≈ 115 with 4 threads

and flat around it: half or twice S* costs about 13% more per sounding. With F ≈ 3,200 s, the best
1-thread run costs ~280 core-s per sounding and the best 4-thread run ~400 (about 100 and 145
core-hours for January 2013's ~1,300 soundings). So, to spend the fewest core-hours:

* **Use 1 thread** (`--threads=1`) and **aim for about 50–120 soundings per run**, taken in time
  order (neighbouring soundings share met and release groups). For January 2013 that is one to
  three consecutive days per run. A 60-sounding run takes about 5 h, a 120-sounding run about 10 h,
  well inside the 24 h limit; above ~150 soundings a 1-thread run gets close to it.
* **Request 16 GB, not 48 GB** (`write_job_script`'s `memory`): every run so far peaked at
  ~10.5 GB, and a 1-core job holding 48 GB blocks other jobs from the node's memory.
* Use 4 threads only when a result is needed quickly (one run in about a quarter of the loop
  time); then aim for 100–200 soundings per run.
* **Combine consecutive days that are small**, **split days with more than ~120 soundings**
  (1 thread) or ~200 (4 threads).
* Limits per run: at most 216 soundings (species particles can carry; not yet reached), 3,000
  field requests, 7 M particles (`MAX_PARTICLES`), 24 h on the `high` queue.

These numbers are good to about ±25%; refine them from the first production months (in particular
a₁ and k₁ in the species layout, measured so far only through the 4-thread runs).
