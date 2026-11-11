# AGENTS.md

Handoff notes for coding agents. Humans: start with README.md.

## Purpose
Rebuild the ACRG NAME satellite-footprint pipeline (run + process) cleanly on JASMIN.
Accuracy is paramount: **reproduce legacy behaviour first, then change it only deliberately.**

## Read before changing anything
- `docs/plan.md`: phases and current status
- `docs/decisions.md`: decisions already made (don't re-litigate without the user)
- `docs/deviations.md`: every intentional difference from legacy output. **Add an entry for any
  change that alters numbers.**
- `dev/notes/`: verified facts about the legacy code, NAME 8.7.1 and the L2 data, with file:line
  references

## Clarity first
The code will be read and checked by scientists, not just run. **Prefer clear over clever or
"perfect" code.**
- Write straightforward, linear code with descriptive names. Avoid deep abstraction, metaprogramming
  and dense one-liners. A little repetition is fine if it reads more easily.
- Every module starts with a docstring that says **where it sits in the workflow** (what comes in,
  what goes out, what step comes next). Every public function's docstring states what it does,
  its inputs and outputs **with units**, and any assumptions.
- Comment the *why*, especially physics, unit conversions and anything that reproduces a legacy
  quirk (cite the legacy file:line and the `docs/deviations.md` ID).
- Keep `docs/workflow.md` (the end-to-end workflow) up to date whenever a step changes.

## Data safety
- `data/` is **read-only** (`chmod -R a-w`). Never make it writable or modify it without the user
  asking. To add new reference data the user (or you, when asked) runs `chmod u+w` on the specific
  folder, adds the files, then restores `a-w`.
- Never point legacy or experimental code at `data/` through symlinks. **Copy** the inputs into a
  new dated scratch folder (`/work/scratch-pw5/$USER/name-satellite/...`) and write outputs only
  there.
- No `rm -rf` on paths built from variables. Use a new folder per attempt; leave cleanup to the user.

## NAME licence
NAME is Met Office software (Crown copyright), used under a licence that restricts its use,
duplication and **disclosure**. This repository is meant to be shareable, so:
- Never commit NAME source code, patches or diffs of it (a diff contains NAME code), executables,
  demo or met-definition files, or met data. The patched build for X20 lives outside the repo.
- Never quote NAME source, or describe its internals with file names, line numbers or internal
  variable names, in code comments or docs. Describe only behaviour a NAME user can observe (input
  syntax, output files, error messages, our measurements): `dev/notes/name-behaviour.md`.
- Notes from reading NAME source go in `~/name-satellite-private/notes/` (outside git).
- NAME files are read where they are installed (e.g. the met declarations on the NAME GWS),
  never copied into the repo.

## Rules
- Tooling: `uv` only (`uv sync`, `uv run ...`, `uv add ...`); never pip-install into the system.
  Python 3.12.
- Before committing: `uv run ruff check . && uv run ruff format --check . && uv run pytest`.
  Pre-commit runs ruff and gitleaks.
- Never commit data, NAME output, `.nc` or `.txt.gz` files, credentials or tokens. Reference data
  live in `data/` (git-ignored); register new ones in `src/name_satellite/refdata.py` with a
  sha256, and in `docs/reference-data.md`.
- Fail loudly: raise on inconsistent input. No silent `return None` (a major legacy problem).
- Match records by explicit IDs (sounding ID, level, source name), never by file or row order.
- Keep functions small and pure where possible; I/O at the edges.
- Compute: small tests and single-day processing are fine on the sci/login node. NAME runs and
  anything heavy go to LOTUS2 via SLURM (see `docs/jasmin.md`). Don't run NAME interactively.
- Commit or push only when the user asks.

## Key paths (JASMIN)
- NAME GWS: `~/name` = `/gws/ssde/j25a/name` (executables and UM met archive)
- Patched NAME 8.7.1 (X20): `/gws/ssde/j25b/acrg/mattrigby/NAME-8.7.1-fateid`, group `gws_name`
  only. Its patch is recorded in a local git repository there; never copy it into this repo.
- ACRG GWS: `/gws/ssde/j25b/acrg`
- Legacy code: `~/code/NAME_scripts`, `~/code/NAME_processing`, `~/code/acrg` (read-only reference)
