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

> **Status: under review, not yet in production.** The legacy workflow is reproduced exactly
> (January 2013 footprints rebuilt bit for bit), and the new workflow runs end to end for a test
> day (1 January 2013). The footprint file format is a draft for the group's review, and the
> domains for production runs are still being decided. See [docs/plan.md](docs/plan.md).

## How it works

**To run it and inspect every file on the way, follow [docs/run-a-day.md](docs/run-a-day.md).**
The full description of each step is in [docs/workflow.md](docs/workflow.md).

1. **Read retrievals** (`observations.py`, `l2_files.py`). Download a day's University of
   Leicester GOSAT or GOSAT-2 CH4 proxy file from CEDA and keep the good soundings in the site's
   box. Each sounding keeps its ID, pressure levels, pressure weights and averaging kernel.
2. **Plan releases** (`prepare.py`, `releases.py`, `release_groups.py`). A short met-only NAME run
   gives the model's surface pressure at each sounding; each sounding's levels are scaled to it.
   Repeat exposures of the same spot share one set of releases. Each group gets one NAME release
   per modelled retrieval level (16 for GOSAT), a layer of air in a box the size of the footprint.
3. **Run NAME** (`name_input.py`, `lotus.py`, `met.py`). Write the NAME input file and run it
   backwards for 30 days on LOTUS2, with UM global met from the NAME group workspace. NAME writes
   a footprint field per release and a record of where each particle left the domain.
4. **Process** (`processing.py`). Turn each release's field into a footprint and each particle
   record into a boundary fate (a side face by height, the top, or still inside). Combine each
   sounding's levels with its own pressure weights and averaging kernel into one column
   footprint, and write a monthly netCDF in the layout of the group's site footprint files
   ([docs/footprint-file-spec.md](docs/footprint-file-spec.md)).

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
uv run pytest                 # run the fast tests
uv run pytest -m ""           # run everything, including slow month-scale checks
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
processing must go to the LOTUS2 batch queue. You need the `gws_name` and `gws_acrg` groups, the
LOTUS2 project `acrg`, and uv in `~/.local/bin` (the job scripts call it there). Working files go
under `/work/scratch-pw5/$USER/name-satellite` unless `NAME_SATELLITE_SCRATCH` /
`NAME_SATELLITE_MET` say otherwise. See [docs/run-a-day.md](docs/run-a-day.md) and
[docs/jasmin.md](docs/jasmin.md).

## Repository layout

```
src/name_satellite/   the Python package
tests/                pytest tests (reference-data tests marked `refdata`)
docs/                 workflow, run-a-day guide, file spec, plan, decisions, deviations, JASMIN
dev/notes/            detailed findings from investigating the legacy code, NAME and data
dev/scripts/          one-off verification scripts used to establish those findings
data/                 local reference data (git-ignored)
```

## Licence and third-party software

This repository is released under the [MIT licence](LICENSE).

- **NAME is not included.** NAME is Met Office software (Crown copyright), available only under
  licence from the Met Office. This repository contains no NAME source code, executables, met
  definitions or met data; it writes NAME input files and reads NAME output, and needs a licensed
  NAME installation (on JASMIN, the NAME group workspace).
- **Satellite data are not included.** The University of Leicester GOSAT and GOSAT-2 proxy XCH4
  products are downloaded from [CEDA](https://data.ceda.ac.uk) when needed, under their own licence
  (Open Government Licence v3).
- **Legacy code.** Parts of the processing reimplement methods from the ACRG legacy code:
  [`acrg`](https://github.com/ACRG-Bristol/acrg) (MIT licence, © Matt Rigby) and `NAME_processing`
  (MIT licence, © hanchawn), and follow the run set-up of `NAME_scripts`.

## Contact

Matt Rigby, University of Bristol (matt.rigby@bristol.ac.uk)
