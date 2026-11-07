# name-satellite

Tools to run Met Office **NAME** back-run footprints for **satellite column retrievals** on
[JASMIN](https://www.jasmin.ac.uk), and to process the NAME output into the footprint files
our inversion systems read.

It replaces the ACRG satellite workflow that was spread across three legacy repositories:

| Legacy code | What it did |
|---|---|
| `acrg.satellite` (in `ACRG-Bristol/acrg`) | Read GOSAT Level-2 files; wrote NAME release CSVs and observation netCDFs |
| `NAME_scripts/satellite` (`ACRG-Bristol/NAME_scripts`) | ksh scripts that built NAME input files and ran NAME |
| `NAME_processing` | Python that turned NAME text output into footprint netCDFs |

> **Status: early development.** Nothing here is ready for production. The current goal is to
> reproduce the legacy workflow exactly, then fix its known problems one at a time.
> See [docs/plan.md](docs/plan.md).

## How it works (target design)

The full step-by-step description is in [docs/workflow.md](docs/workflow.md).

1. **Read retrievals.** Read satellite Level-2 files (first: University of Leicester GOSAT and
   GOSAT-2 CH4 proxy) into a standard table of soundings, each with a unique ID, its pressure
   levels, pressure weights and averaging kernel.
2. **Plan releases.** For each sounding, define a NAME particle release for each retrieval
   level up to a chosen pressure ceiling.
3. **Run NAME.** Write NAME input files and run them on the LOTUS2 batch cluster, using UM
   global met from the NAME group workspace.
4. **Process.** Read NAME's fields and particle files into footprints (time-integrated or
   hourly), boundary-crossing particle histograms and release-point met. Combine the per-level
   footprints into one column footprint using the pressure weights and averaging kernel. Write
   the netCDF format the inversion code expects.

## Code style

The code is meant to be read and checked by scientists. We aim for **clarity over perfect code**:
plain, well-named functions, and docstrings that say where each step sits in the workflow, with
units.

## Getting started

You need [uv](https://docs.astral.sh/uv/). On JASMIN it can be installed into `~/.local/bin`.

```bash
git clone <repo-url> name-satellite
cd name-satellite
uv sync                       # creates .venv with Python 3.12 and all dependencies
uv run pre-commit install     # ruff + gitleaks checks on every commit
uv run pytest                 # run the tests
```

Useful commands:

```bash
uv run ruff check .           # lint
uv run ruff format .          # format
uv run pre-commit run --all-files
```

### Reference data

Tests that compare against legacy output need reference files that are **not** in git (they are
large and some are not ours to publish). Put them in `data/`, or set
`NAME_SATELLITE_REFDATA=/path/to/them`. Tests skip automatically when the files are missing.
[docs/reference-data.md](docs/reference-data.md) lists the files, where they came from and their
checksums. `data/` is kept **read-only** (`chmod -R a-w data`) so nothing can change it by accident.
To add files, make just that folder writable (`chmod u+w`), add them, then make it read-only again.

## Running on JASMIN

Small tests and single-day processing can run on a JASMIN sci node. NAME runs and month-scale
processing must go to the LOTUS2 batch queue. See [docs/jasmin.md](docs/jasmin.md).

## Repository layout

```
src/name_satellite/   the Python package
tests/                pytest tests (reference-data tests marked `refdata`)
docs/                 plan, decisions, deviations from legacy, reference data, JASMIN notes
dev/notes/            detailed findings from investigating the legacy code, NAME and data
dev/scripts/          one-off verification scripts used to establish those findings
data/                 local reference data (git-ignored)
```

## Contact

Matt Rigby, University of Bristol (matt.rigby@bristol.ac.uk)
