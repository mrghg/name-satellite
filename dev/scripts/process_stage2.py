"""Process and compare the stage 2 NAME 8.7.1 validation runs (1-7 Jan 2013, GOSAT-BRAZIL).

Runs (prepare_stage2.py), all processed with exactly the legacy processing
(``legacy_pipeline.process_legacy_inputs``); release-point met from the legacy 2018 run, which has
all days (same release points and times):

    day_<YYYYMMDD>  2-7 Jan, one run per day, legacy window   ->  compared with production
    batch_0107      1-7 Jan in one run, corrected window      ->  split into days (batching.py),
                                                                  compared with the daily runs

Yardstick: the stage 1 noise test (run B vs B2, which differ only in their random seed), read
from the stage 1 folder if it has been processed. Differences no larger than that are noise.

Batch vs daily: the daily runs use the legacy window, the batched run the corrected one. They
differ for soundings released later than the day's first (X18), so "identical?" is checked only
for the day's EARLY soundings (released within 10 min of the first). For 1 Jan the batch is
compared with stage 1 run A, which has the corrected window too.

    uv run python dev/scripts/process_stage2.py \
        /work/scratch-pw5/mattrigby/name-satellite/stage2_20261109 \
        /work/scratch-pw5/mattrigby/name-satellite/stage1_20261109
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from name_satellite.batching import split_batched_output
from name_satellite.legacy import read_legacy_met_folder
from name_satellite.legacy_pipeline import (
    legacy_particles_from_fate,
    process_legacy_inputs,
    write_legacy_fp_file,
)
from name_satellite.name_output import read_fields_folder, read_name_fate_file
from name_satellite.refdata import ref_path
from name_satellite.validation import histogram_difference, shape_difference, total_ratio

LEGACY_RUN = ref_path("conc2018_input").parents[1]
RUN_NAME = "SOUTHAMERICA_GOSAT-BRAZIL_column"
DAYS = [f"201301{d:02d}" for d in range(1, 8)]
EDGES = "nesw"
FACTORS = [1, 4, 16]
EARLY_MINUTES = 10


def finished(run):
    """True if the NAME run wrote its fields file(s) and its error file is empty."""
    error = run / "inputError.txt"
    fields = list((run / "output").glob("Fields_grid1_C1_T1_*.txt"))
    return len(fields) == 1 and error.exists() and error.stat().st_size == 0


def read_run(run):
    # A batched run may write several fields files (name_input.fields_output_group).
    fields = read_fields_folder(run / "output")
    particles = legacy_particles_from_fate(
        *read_name_fate_file(run / "output" / "Particles_C1.txt")
    )
    return fields, particles


def process_day(fields, particles, day, path):
    """Legacy processing of one day; writes the fp file and returns it as written."""
    met = read_legacy_met_folder(LEGACY_RUN / "Met" / f"Met_{RUN_NAME}_{day}")
    obs_file = LEGACY_RUN / "Observations" / f"gosat-fts_gosat_{day}_ch4-column.nc"
    obs = xr.open_dataset(obs_file).load()
    result = process_legacy_inputs(fields, particles, met, obs, release_seconds=60.0)
    write_legacy_fp_file(result, path)
    return xr.open_dataset(path).load()


def arrays(ds):
    """Footprints (sounding, lat, lon) and edge histograms {edge: (sounding, position, height)}."""
    fp = ds["fp"].transpose("time", "lat", "lon").values
    hist = {}
    for e in EDGES:
        along = "lon" if e in "ns" else "lat"
        hist[e] = ds[f"particle_locations_{e}"].transpose("time", along, "height").values
    return fp, hist


def summary(label, x, ref):
    """One row of medians comparing processed footprints x against ref (same soundings)."""
    fx, hx = arrays(x)
    fr, hr = arrays(ref)
    row = {
        "pair": label,
        "soundings": len(fx),
        "total ratio median": np.median(total_ratio(fx, fr)),
        "total ratio std(log)": np.log(total_ratio(fx, fr)).std(),
    }
    for factor in FACTORS:
        row[f"shape {factor}x"] = np.median(shape_difference(fx, fr, factor))
    row["hist"] = np.median(histogram_difference(hx, hr))
    return row


def early_soundings(ds):
    """Indices of soundings released within EARLY_MINUTES of the day's first."""
    times = pd.to_datetime(ds["time"].values)
    return np.flatnonzero(times - times.min() <= pd.Timedelta(minutes=EARLY_MINUTES))


def main(work, stage1):
    work, stage1 = Path(work), Path(stage1)
    production = xr.open_dataset(ref_path("production_201301"))
    rows = []

    # Yardstick from stage 1: B vs B2 (seed only).
    b_file = stage1 / "run_B" / "fp_20130101_legacy_processing.nc"
    b2_file = stage1 / "run_B2" / "fp_20130101_legacy_processing.nc"
    if b_file.exists() and b2_file.exists():
        rows.append(
            summary("NOISE 0101 B vs B2", xr.open_dataset(b_file), xr.open_dataset(b2_file))
        )
    else:
        print("Stage 1 noise test not processed yet: no yardstick row")

    # Daily runs vs production.
    daily = {}
    for day in DAYS[1:]:
        run = work / f"day_{day}"
        if not finished(run):
            print(f"{run.name}: not finished, skipped")
            continue
        fields, particles = read_run(run)
        daily[day] = process_day(fields, particles, day, run / f"fp_{day}_legacy_processing.nc")
        prod_day = production.sel(time=daily[day]["time"].values)
        rows.append(summary(f"{day} daily vs production", daily[day], prod_day))

    # Batched run, split into days.
    run = work / "batch_0107"
    if finished(run):
        plan = pd.read_csv(run / "release_plan.csv", dtype={"day": str})
        fields, particles = read_run(run)
        run_a = stage1 / "run_A" / "fp_20130101_legacy_processing.nc"
        for day in DAYS:
            day_fields, day_particles = split_batched_output(fields, particles, plan, day)
            batch = process_day(
                day_fields, day_particles, day, run / f"fp_{day}_legacy_processing.nc"
            )
            if day == DAYS[0]:
                reference, name = (xr.open_dataset(run_a).load() if run_a.exists() else None), "A"
            else:
                reference, name = daily.get(day), "daily"
            prod_day = production.sel(time=batch["time"].values)
            rows.append(summary(f"{day} batch vs production", batch, prod_day))
            if reference is None:
                print(f"{day}: no {name} run to compare the batch with")
                continue
            rows.append(summary(f"{day} batch vs {name}", batch, reference))
            early = early_soundings(batch) if name == "daily" else np.arange(batch.sizes["time"])
            same = np.array_equal(
                batch["fp"].isel(time=early).values, reference["fp"].isel(time=early).values
            )
            print(
                f"{day}: batch fp identical to {name} for {len(early)} compared soundings: {same}"
            )
    else:
        print(f"{run.name}: not finished, skipped")

    table = pd.DataFrame(rows).set_index("pair")
    pd.set_option("display.width", 220)
    print(table.round(3).to_string())
    table.to_csv(work / "stage2_comparison.csv")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
