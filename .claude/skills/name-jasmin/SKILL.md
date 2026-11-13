---
name: name-jasmin
description: Practical knowledge for running the Met Office NAME III dispersion model (8.7.1) on JASMIN and LOTUS2, learned from building the ACRG satellite-footprint pipeline - input-file syntax and its traps, UM global met versions and the archive, Slurm job set-up, output and particle-fate files and their quirks, memory and run-time measurements. Use this skill whenever the task involves NAME runs, NAME input files, NAME output (fields files, Particles_C1.txt, inputLog/inputError/inputTiming), UM met (UMG_Mk6..Mk12, met restore), back-runs or footprints made with NAME, or submitting NAME jobs to LOTUS2, even if the user does not name the skill or says only "the model" or "the back-runs".
---

# NAME III on JASMIN

What has been learned, by running it, about the Met Office NAME III model (version 8.7.1) on
JASMIN's LOTUS2 cluster. Everything here is **observable behaviour**: input syntax, output files,
error messages and our own measurements.

## The licence comes first

NAME is Crown-copyright Met Office software used under a licence that restricts use, duplication
and disclosure. Anything you write that others may see (a git repository, docs, slides, issues)
must therefore:

- never contain NAME source code, patches or diffs of it, executables, demo or met-definition
  files, or met data;
- never quote NAME source or describe its internals (source file names, line numbers, internal
  routine or variable names). Describe only what a NAME user can observe: input syntax, output
  files, error messages, measurements.

If you read NAME source to understand something, keep those notes private (for the ACRG
satellite work: `~/name-satellite-private/notes/`, mode 700, outside git) and translate the
finding into observable behaviour before writing it anywhere shared. Read NAME files where they
are installed; do not copy them into a repository.

## Before writing or submitting a run

1. **Test the layout small first.** Some failures only appear at the *end* of a run (e.g. an
   output line too long, after 8 hours). Run the same input with `Max Travel Time` of 1 h on the
   `debug` queue first: minutes, and it checks syntax, grids, output groups and memory.
2. **Restore the met before any NAME job runs, in one job.** Never let several NAME jobs restore
   met into the same folder at once: they read half-unzipped files and lose particles.
3. **Judge a run by `<input>Error.txt`, not the exit code.** NAME exits 0 even when it lost
   particles. Make job scripts fail when the error file is not empty, so dependent jobs stop.
4. **Threads must match three ways:** Slurm `--cpus-per-task`, `OMP_NUM_THREADS` and the input's
   `Threads`. NAME will not use more threads than OpenMP offers. 4 threads ≈ 80% efficient.
5. **Multi-CPU jobs need QoS `high`** (or `debug`, up to 8 CPUs, 1 h). `short`, `standard` and
   `long` allow one CPU. Charge to the right account (`--account=acrg` for ACRG work).
6. **Don't poll Slurm in tight loops.** Check `squeue`/`sacct` now and then; if you must wait in
   a loop, no more than once a minute (JASMIN asks users not to poll Slurm in tight loops).
7. **`Random Seed = Fixed`** makes runs repeatable; `Random` gives a new sequence each run (no
   numeric seed is documented). The docs warn threaded runs may differ slightly; our 1/4/8-thread
   runs were bit-identical. Measure particle noise as Fixed vs Random (or two Random runs).
8. **Check the NAME documentation** (below) before relying on memory for syntax: it is
   authoritative for input and output formats; our notes add what it doesn't cover.

## Where to look next

Read only the file the task needs:

| Task | Read |
|---|---|
| Writing or debugging an input file (blocks, domains, grids, sources, requests, limits) | `references/input-syntax.md` |
| Choosing or restoring UM met; dates, versions, changeovers, archive coverage | `references/met.md` |
| Reading output: fields files, particle-fate records and their duplicates, point output | `references/outputs.md` |
| Job scripts, environment, queues, storage, measuring memory | `references/jasmin-lotus.md` |
| How long a run takes, how much memory, what drives cost | `references/performance.md` |

## NAME's own documentation

`/gws/ssde/j25a/name/docs/namedoc_v8_6/docs/` (PDFs, for NAME 8.6; we run 8.7.1, so Mk12 and some
behaviour are newer than the docs). Read with `pdftotext -layout <file> <scratch>.txt` and grep.
They are Crown copyright and "not released for publication": paraphrase and cite (document +
section), never copy passages into anything shared.

| Task | Where |
|---|---|
| Input syntax, every block and column | B01 Input (§6–7 and the Appendix on syntax) |
| Grids, domains, back-runs, sources | A01 User Guide §1.2.6, §1.2.7.3, §1.5.4 |
| Field requests, averaging, output formats | B01 "Output Requirements – Fields"; A01 §1.4.1, §1.4.4.1 (Tables 1.12–1.13); B02 §2 |
| Particle and fate output | B01 "Sets of Particle/Puff Information"; A01 §1.4.2, §1.4.4.2, §1.2.6.5 |
| Met set-up, UM versions | A03 Using NWP data §2.1, Table 8, App. A.1; A01 §1.6.2; B01 NWP Met / Flow Module Instances |
| Threads, run time | A01 §1.2.7.7, §1.8.4, §1.8.5.2 |
| Error messages, message levels | B03 (old but still relevant); A01 §1.2.7.8 (Message Controls) |

## Key paths (JASMIN)

| What | Path |
|---|---|
| NAME group workspace (executables, UM met archive, demo scripts) | `/gws/ssde/j25a/name` (`~/name`; needs Unix group `gws_name`) |
| NAME 8.7.1 | `/gws/ssde/j25a/name/code/NAME-8.7.1/Executables_Linux/nameiii_64bit_par.exe` |
| ACRG-patched 8.7.1 (fate records carry particle index and fate code) | `/gws/ssde/j25b/acrg/mattrigby/NAME-8.7.1-fateid` (group `gws_name` only) |
| NAME team's 8.7.1 demo scripts and met declarations | `/gws/ssde/j25a/name/examples/DemoScripts_v8_7_1/` |
| UM global met archive | `/gws/ssde/j25a/name/met_archive/Global/UMG_Mk*PT/` |
| UM topography | `/gws/ssde/j25a/name/data/UMTopogData` |
| Scratch for runs and restored met | `/work/scratch-pw5/$USER/...` (files deleted 28 days after last access) |

The ACRG satellite pipeline that these notes come from is `~/code/name-satellite`
(github.com/mrghg/name-satellite): `dev/notes/name-behaviour.md` and `docs/jasmin.md` there hold
the same knowledge with the evidence (job numbers, dates), and its `name_input.py`, `lotus.py`
and `met.py` are working examples of input files, job scripts and met handling.

## Keeping this skill current

When you learn something new about NAME's behaviour (a limit, an error message, a measurement),
add it to the matching reference file with the date and how it was found, in observable terms.
