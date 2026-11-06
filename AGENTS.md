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
- ACRG GWS: `/gws/ssde/j25b/acrg`
- Legacy code: `~/code/NAME_scripts`, `~/code/NAME_processing`, `~/code/acrg` (read-only reference)
