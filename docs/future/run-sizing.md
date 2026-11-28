# Sizing future-method runs (draft, 28 Nov 2026)

How many soundings to put in one NAME run, for the fewest core-hours. Built from the January and
August 2013 test runs (10-day back-runs, mostly 4 threads, `future-method` set-up); good to about ±25%,
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

## What a run occupies on a node

A LOTUS2 `standard` node has 192 cores and 1.54 TB: about **8 GB per core**. A job takes its
cores *or* its memory, whichever is the larger share of the node, so count its cost as
core-equivalents = max(threads, memory requested / 8 GB), times its run time. Our Mk6/Mk7 runs
peak at ~10.5 GB (16 GB requested):

| Threads | Core-equivalents | 7 Jan part 1 (s) | Cost (core-equivalent-s) |
|---|---|---|---|
| 1 | 2 (memory) | 12,055 | 24,100 |
| **2** | **2** | 7,488 | **15,000** |
| 4 | 4 | 5,000–6,800 | 20,000–27,000 |

One thread looks cheapest in cores alone, but a node full of 1-thread, 16 GB jobs runs out of
memory with half its cores idle. **Two threads match the node's ratio of memory to cores**, so
nothing is stranded whatever gap the job fills, and they are the default (`THREADS = 2`). (If
LOTUS2's fair-share counts cores only, 1 thread would look cheaper on the account while taking more
of the cluster.)

## Best run size (species layout, 2 threads)

With S soundings in a run, p particles per sounding (January 2013: p ≈ 12,800) and c
core-equivalents, the cost per sounding is c·F/S + c·a_c·p + 3·c·k_c·p·S, where a_c and k_c are the
values with that many threads (with 2 threads, about 1.9× faster than 1: c·a_c ≈ 14,700 s and
c·k_c ≈ 23 s per million particles, with F ≈ 3,200 s). It is smallest at

    S* = √(c·F / (3·c·k_c·p))   ≈ 85 soundings with 2 threads (≈ 120 with 4)

and flat around it (half or twice S* costs about 13% more). The best 2-thread run costs ~340
core-equivalent-s per sounding, against ~410 with 4 threads and ~560 with 1: about 120
core-equivalent-hours for January 2013's ~1,300 soundings. So:

* **Aim for about 50–150 soundings per run (best ~85)**, taken in time order (neighbouring soundings
  share met and release groups). For January 2013 that is two consecutive days per run. An
  85-sounding run takes about 4–5 h.
* **Combine consecutive days that are small; split days with more than ~150 soundings.**
* Use 4 threads only when a result is needed quickly.
* Limits per run: at most 216 soundings (species particles can carry; not yet reached), 3,000
  field requests, 7 M particles (`MAX_PARTICLES`), 24 h on the `high` queue.

These numbers are good to about ±25%; the 2-thread figure is one run. Refine them from the first
production months (in particular in the species layout, so far measured only with 4 threads).

## Memory and met version (outstanding: not yet measured beyond Mk7)

Peak memory is mostly met (an unbounded run holds global met: +5.6 GB over a bounded one, see the
`name-jasmin` skill), plus the output grids. Mk6 and Mk7 have the same grid, so the same memory.
Later versions are finer (first file of each version, 28 Nov 2026):

| Version | Grid (whole globe) | Points vs Mk6 | One met time, as stored |
|---|---|---|---|
| Mk6, Mk7 | 1024 × 769 (~25 km) | 1× | 1.6 GB |
| Mk8, Mk9 | 1536 × 1152 (~17 km) | 2.2× | 3.6 GB |
| Mk10–Mk12 | 2560 × 1920 (~10 km) | 6.0× | 2.4 GB (packed) |

NAME holds the fields unpacked, so memory should follow the grid points, not the file sizes. If
the ~6–9 GB of met in today's ~10.5 GB scales with the points: **Mk8/9 roughly 17–25 GB, Mk10–12
roughly 40–60 GB**. Met reading (F) will also grow, by 1.5–2× in bytes read and more if unpacking
costs time. Both change the sizing: with ~20 GB the node ratio is met by 2–3 threads, with ~50 GB by
about 6, and a larger F favours larger runs. Until measured, jobs on these versions request 48 GB
(`JOB_MEMORY`); a 1-hour test run on each grid would measure the memory (it is reached at the
start), and a full run the read time.
