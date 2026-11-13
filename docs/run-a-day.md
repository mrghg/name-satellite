# Running one day end to end, and checking every file

This guide runs the whole new-run workflow for one day and shows what to look at after each step.
The example is the one we tested: **GOSAT-1 v10.0 over GOSAT-BRAZIL, 1 January 2013** (46
soundings, 22 release groups, 352 NAME releases). [workflow.md](workflow.md) explains what each
step does and why; the module docstrings give the details.

Times below are from our test (November 2026). Expect about **1.5 hours** of LOTUS2 time for the
NAME run, a few minutes for everything else, and up to a few hours for the met restore if you
don't reuse ours (step 2).

## 0. Before you start

You need:

- a JASMIN account in the Unix groups **`gws_name`** (the NAME group workspace: NAME itself, the
  UM met archive and our patched NAME build) and **`gws_acrg`**; check with `groups`;
- the LOTUS2 project **`acrg`**, which the jobs are charged to (`--account=acrg`);
- [uv](https://docs.astral.sh/uv/) installed in **`~/.local/bin/uv`** (the job scripts call it
  there, because batch jobs don't load your shell's `PATH`).

Then, on a sci node (e.g. `sci-vm-01`):

```bash
git clone git@github.com:mrghg/name-satellite.git ~/code/name-satellite   # private: ask Matt for access
cd ~/code/name-satellite
uv sync
uv run pytest -q          # ~25 s; tests that need reference data skip without them
```

Working files go in scratch, not your home folder (a day of NAME output is about 0.6 GB; scratch
files are deleted 28 days after their last access). By default everything goes under
`/work/scratch-pw5/$USER/name-satellite`. Two environment variables change that:

| Variable | Default | What it is |
|---|---|---|
| `NAME_SATELLITE_SCRATCH` | `/work/scratch-pw5/$USER/name-satellite` | L2 download cache, restored met (unless the next variable is set) |
| `NAME_SATELLITE_MET` | `$NAME_SATELLITE_SCRATCH/met` | where the UM met is restored, one subfolder per met version |

The met for a 30-day back-run of January 2013 is 553 GB once unzipped. Ours is already restored
and readable by everyone, so reuse it rather than make a second copy:

```bash
export NAME_SATELLITE_MET=/work/scratch-pw5/mattrigby/name-satellite/met
```

Set these in the same shell you submit from: Slurm jobs inherit the environment they are submitted
with.

## 1. Prepare the day (sci node, ~1 minute)

```bash
WORK=/work/scratch-pw5/$USER/name-satellite/review_$(date +%Y%m%d)
uv run python -m name_satellite.prepare "GOSAT-1 v10.0" GOSAT-BRAZIL 2013-01 $WORK --days=2013-01-01
```

This downloads the day's L2 file from CEDA (into the L2 cache), keeps the good soundings in the
site's box, and writes the jobs. It refuses to write into a month folder that already exists, so
use a new `WORK` for each attempt. Add `--test` for a quick check of the set-up first: particles then
travel for 1 hour only, on the debug queue (minutes instead of 1.5 hours).

The month folder is `$WORK/GOSAT-BRAZIL_GOSAT-1_v10.0_201301/`. Check:

| File | What to look at |
|---|---|
| `summary.csv` | one row per day: L2 soundings in the file (2,179), soundings kept (46), NAME releases (352), release groups (22), met version and the met period needed |
| `20130101/obs_GOSAT-BRAZIL_20130101.nc` | **the processed L2 observations**, one entry per sounding (below) |
| `20130101/day.json` | product, site, day, threads, test flag (what was asked) |
| `20130101/surface_pressure/soundings.csv`, `input.txt` | the met-only NAME run: one sample point per sounding, at its position and minute |
| `restore_met.sbatch`, `20130101/surface_pressure/job.sbatch`, `20130101/job.sbatch`, `submit.sh` | the three jobs and the script that chains them |

The observation file (`ncdump -h`, or in Python):

```python
import xarray as xr

obs = xr.open_dataset(
    f"{WORK}/GOSAT-BRAZIL_GOSAT-1_v10.0_201301/20130101/obs_GOSAT-BRAZIL_20130101.nc"
)
obs  # dims sounding_id (46) and lev (20, surface first)
obs.attrs  # product, version, source L2 file, bias correction, time resolution
obs.pressure_levels[0]  # hPa; pressure_weights and xch4_averaging_kernel on the same levels
```

`sounding_id` is the L2 exposure ID (D18). Times are whole minutes, rounded down, because that is
all the v10 files store (X22). XCH4 has the GOSAT-1 bias correction applied (see the attribute).

## 2. Submit (sci node), then wait

```bash
bash $WORK/GOSAT-BRAZIL_GOSAT-1_v10.0_201301/submit.sh
squeue -u $USER          # check now and then; please don't poll in a loop
```

`submit.sh` submits three jobs, each starting only if the one before succeeded (`afterok`):

1. **met restore** (`restore-<job>.out` in the month folder): makes every met file of the period
   available in the met folder. Seconds if they are already there; hours for a fresh month of
   Mk6 met.
2. **surface-pressure run** (`20130101/surface_pressure/`, stock NAME 8.7.1, seconds to minutes):
   NAME's surface pressure at each sounding, plus the pressure at the 20 km domain top.
3. **main job** (`20130101/`, patched NAME 8.7.1 on 4 threads, ~82 minutes here): first writes
   the NAME input, then runs NAME.

A job script exits non-zero if NAME fails *or* writes anything to its error file
(`inputError.txt`), so a failed day stops its later jobs. When a job fails, read its
`slurm-<job>.err`/`.out` and then NAME's `inputError.txt` and the end of `inputLog.txt`.

## 3. Check the NAME input (written by the main job before NAME starts)

In `20130101/`:

| File | Rows | What to check |
|---|---|---|
| `surface_pressure.csv` | one per sounding | `retrieval_hpa`, `model_hpa`, their `ratio` (the job stops outside 0.9–1.1; every level is scaled by it, D19) and `domain_top_hpa`, the model pressure at the 20 km top: the 16 modelled layers must lie below it, i.e. at higher pressure (D21) |
| `release_groups.csv` | one per sounding | `sounding_id` → `release_group` (the group's earliest sounding ID). Here 12 groups of 3 and 10 single soundings |
| `release_plan.csv` | one per NAME release | `source_name` (`SourceID_<n>`), `release_group` (the group's ID, i.e. its earliest sounding's ID; `point_number` is the group's number), `level_number` 1–16, release window (60 s; `release_end` is earlier because NAME runs backwards), centre `lon`/`lat` and box `dlon`/`dlat` (degrees, full widths), layer centre `z_pa` and thickness `dz_pa` (Pa), `particles_per_second` |
| `input.txt` | — | the NAME input file (below) |

Quick look, e.g. that each group has levels 1–16 and the layers are stacked without gaps:

```python
import pandas as pd

plan = pd.read_csv("release_plan.csv")
plan.groupby("release_group").level_number.agg(["min", "max", "count"])
g = plan[plan.release_group == plan.release_group[0]].sort_values("level_number")
(g.z_pa - g.dz_pa / 2).values[:-1] - (g.z_pa + g.dz_pa / 2).values[1:]  # ~0: layers touch
```

`input.txt` is in NAME's usual block format: a heading line ending in `:`, a header row, then one
row per item. The blocks to check are `Domains:` (the dispersion domain: lon/lat bounds, 20,000 m
agl top, 720 h maximum travel), `Sources:` (one row per release, matching `release_plan.csv`),
`Output Requirements - Fields:` (one time-integrated mixing-ratio field per source, on the output
grid), `Output Requirements - Sets of Particle/Puff Information:` (the particle-fate records, with
heights in m agl), and, at the end, the release-point met requests (one `ReleaseMet_S<nnnn>`
group per release group). List the blocks with `grep -n ":$" input.txt`.

## 4. Check the NAME output

In `20130101/`: `inputLog.txt` should end with "Run completed successfully", `inputError.txt`
should be empty, and `inputTiming.txt` gives the run time. In `20130101/output/` (~600 MB):

| File | What it is |
|---|---|
| `Fields_grid1_C1_T1_<start>.txt` (and `Fields_grid1_2_...` when there are more than 2,000 releases) | one column per release: the mixing ratio at the release integrated over the 30 days (`ppm s`) on the output grid (190 × 357 cells of 0.352° × 0.234°). The header names the source of each column (row starting `SourceID_`), the grid origin and size; data rows are `X Index, Y Index, X, Y, value...`, and all-zero rows are left out |
| `Particles_C1.txt` | one record per particle when it stops: source number, travel time (h), position (lon, lat, height m agl), particle index within its source, fate code (8 = left the domain, 16 = still inside at 720 h). About 545,000 records here |
| `ReleaseMet_S<nnnn>_C1.txt` | met at each release group's lowest release, at its time: temperature, pressure, wind, boundary-layer depth, surface pressure |

These are big text files: use `head -40` (fields) or `sed -n 11,13p` (particles), not an editor.
Two known features of the particle file are handled in processing and described in
`dev/notes/name-behaviour.md`: about 1.7% of particles appear twice when they leave (X20; the
first record is kept), and a few exits sit exactly on the domain edge as printed.

## 5. Process into the footprint file (sci node, ~20 s, 1.5 GB)

```bash
uv run python -m name_satellite.processing $WORK/GOSAT-BRAZIL_GOSAT-1_v10.0_201301 $WORK/processed --days=20130101
```

Without `--days` it processes every day folder of the month (do that on LOTUS2 for a full month).
It stops at a day whose NAME run failed, so check `sacct` or the day's `inputLog.txt` first.
It stops with an error if anything doesn't add up: sources, groups and soundings must match by ID,
every release's particle fates must sum to 1, and each sounding's to its column weight.

The output is `$WORK/processed/GOSAT1v10-BRAZIL-column_NAME_UMG_SOUTHAMERICA_inert_201301.nc`
(2.7 MB here). [footprint-file-spec.md](footprint-file-spec.md) defines every variable. To look
at it:

```python
import numpy as np, xarray as xr

fp = xr.open_dataset(
    f"{WORK}/processed/GOSAT1v10-BRAZIL-column_NAME_UMG_SOUTHAMERICA_inert_201301.nc"
)
fp.attrs  # provenance: model, met, deviations and decisions applied, commit
srr = fp.srr.isel(footprint=fp.footprint_index)  # each sounding's footprint, (mol/mol)/(mol/m2/s)
srr.isel(sounding_id=0).plot(robust=True)  # a map; most sensitivity is near the sounding

# Every particle is counted once: faces + top + inside = the column weight, per sounding.
w = fp.column_weight.sum("release_level")
fates = fp.particle_fraction_inside + fp.particle_locations_top.sum(("latitude", "longitude"))
for f in "nesw":
    fates = fates + fp[f"particle_locations_{f}"].sum(
        ("height", "longitude" if f in "ns" else "latitude")
    )
np.abs(fates.values[fp.footprint_index.values] / w.values - 1).max()  # < 1e-4 (rounding)
```

What we saw for this day: fates/column weight within 3e-5 of 1; median shares of the column
leaving west 66%, south 9%, east 7%, north 0.3%, top 0.1%, still inside 1.2%; the shares by face
match the legacy production footprints of the same day to 0.003 (plan.md, Phase 6).

## Where to read more

- [workflow.md](workflow.md): every step, the modules that do it, and the decisions behind it.
- [decisions.md](decisions.md) (D-numbers) and [deviations.md](deviations.md) (X-numbers): what
  differs from the legacy pipeline, and why.
- [jasmin.md](jasmin.md): LOTUS2 queues, storage, paths.
