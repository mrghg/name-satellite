# UM global met for NAME on JASMIN

Checked 13 Nov 2026 by listing `/gws/ssde/j25a/name/met_archive/Global/`. Working code:
`~/code/name-satellite/src/name_satellite/met.py` (`MET_VERSIONS`, `met_version_for_period`,
`required_met_files`, `missing_met_files`, `restore_met_files`).

## Versions (59-level, partitioned "PT" files: 14 parts, 3-hourly)

| Version | Official use (NAME team's 8.7.1 demo dates) | Archive complete from | In archive as |
|---|---|---|---|
| UMG_Mk6PT | 2010-03-09 – 2013-04-29 | (official start) | gzipped |
| UMG_Mk7PT | 2013-04-30 – 2014-07-14 | 2013-04-11 | gzipped |
| UMG_Mk8PT | 2014-07-15 – 2015-08-24 | 2014-06-13 | gzipped |
| UMG_Mk9PT | 2015-08-25 – 2017-07-10 | 2015-07-21 | gzipped |
| UMG_Mk10PT | 2017-07-11 – 2022-05-03 | 2017-06-02 | plain `.pp` |
| UMG_Mk11PT | 2022-05-04 – 2026-01-20 | 2022-04-02 | plain `.pp` |
| UMG_Mk12PT | 2026-01-21 – current | 2025-12-02 | plain `.pp` (8.7.1 only) |

- Official dates come from the NAME team's `BackRuns_OneSite_ByDay.sh` (8.7.1 demo scripts).
- Every file NAME would ask for in each official period is present (checked ~1.2 M files).
- From Mk7 there are two files per time and part (`_I_` and `_M_`). File names come from the met
  definition file, so build the list from it rather than guessing.
- Older global met also exists, not set up in our code: UMG_Mk5 (52 levels, not partitioned),
  GLOUM6pp, GLOUM6, GLOUM5, GLOH. NAME 8.7.1 has
  definitions for Mk5 and UM6G. Official dates (A03 Table 8): UM6G pp 01/2009 – 10 Nov 2009,
  Mk5 10 Nov 2009 – 9 Mar 2010; the archive holds Mk5 files from 6 Oct 2009. GOSAT starts April
  2009, so its first year needs these.
- Mk6 also has 52-level partition definitions (L52PT) that ignore levels above ~20 km, for speed
  (A03 App. A.1). Using them would change results: log a deviation first.
- Mk6 partitions are 1024 × 45 points (~25 km) each (A03 App. A.1).
- Met top ≈ 29 km.

## A run uses one version
A back-run must use one met version throughout. Near a changeover, use the version in use on the
run's **last** day of met for the whole run: each newer version is archived weeks before its
official start (Mk8 32 days, Mk9 35, Mk10 39, Mk11 32, Mk12 50), enough for a 30-day back-run.
Exception: Mk7 is archived only 19 days early, so 30-day back-runs for days 29 Apr – 11 May 2013
cannot use one version (shorten them, or switch version during the run). The legacy ACRG scripts
instead shortened back-runs by hand near changeovers and stopped at Mk10.

## Restoring met
- NAME's on-demand restore (the `MetRestore_JASMIN.sh` script named in the declarations) unzips
  gzipped versions into the met folder and links plain ones.
- Sizes: Mk6 one day ≈ 3 GB gzipped, ≈ 12 GB unzipped; a 30-day back-run of one month ≈ 390–550 GB.
  Mk10–Mk12 cost no space (links).
- **Restore everything a batch of runs needs in one job before any NAME run starts.** Several
  NAME jobs restoring into one folder at once read half-written files: "ERROR: in reading met
  data ... Particles and/or puffs are being lost" in `<input>Error.txt`, exit code still 0.
  Safe pattern: decompress to `<name>.partial`, then rename; keep files already present only if
  their size matches.
- One shared met folder per version in scratch serves all runs; set `Delete Met? = No` (Yes
  deletes working copies after reading). Never point the met folder at the archive itself
  (A01 §1.6.2).
- `Update on Demand? = Yes` must be set in both the met module instance and the flow module
  instance; NAME then reads only the partitions it needs (A01 §1.6.2).
- B03: with global met, lost-particle errors most likely mean missing met files.
- Scratch deletes files 28 days after last access.

## Memory
The declarations use each met file's whole global domain. Peak memory is dominated by met and
output grids, not particles: see performance.md (bounded domain ~4–5 GB; unbounded +5.6 GB).
