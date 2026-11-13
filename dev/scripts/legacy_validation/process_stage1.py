"""Process and compare the stage 1 NAME 8.7.1 validation runs (1 Jan 2013, GOSAT-BRAZIL).

The runs (prepare_stage1.py, prepare_stage1_noise.py) go through exactly the legacy processing
(``legacy_pipeline.process_legacy_inputs``), with two stand-ins for inputs the new runs don't make:

* release-point met: the legacy 2019 met-extraction files (same release points and times);
* particles: the NAME 8.7.1 fate file (``legacy_particles_from_fate``).

Writes ``<work>/run_<X>/fp_20130101_legacy_processing.nc`` and prints:

1. Per sounding:
       A vs B           effect of the integration window (X18)
       B vs production  NAME 8.7.1 vs NAME 7.2 with the same window
2. If run_B2 has finished, the noise test. B2 differs from B only in its random seed, so B - B2 is
   pure particle noise. If NAME 8.7.1 reproduces legacy, B - production should be no larger:
       * scatter of the per-sounding total ratios;
       * footprint shape differences at several grid coarsenings (validation.shape_difference);
       * edge histogram differences;
       * per-release footprint totals by level (from the raw fields), including the mean
         B/production ratio with its standard error (a bias test).

    uv run python dev/scripts/legacy_validation/process_stage1.py /work/scratch-pw5/mattrigby/name-satellite/stage1_20261109
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from name_satellite.name_output import read_fields, read_name_fate_file
from name_satellite.validation.compare import (
    histogram_difference,
    near_field_ratio,
    shape_difference,
    total_ratio,
)
from name_satellite.validation.legacy_files import read_legacy_met_folder, read_name_input_sources
from name_satellite.validation.legacy_pipeline import (
    legacy_particles_from_fate,
    process_legacy_inputs,
    write_legacy_fp_file,
)
from name_satellite.validation.refdata import ref_path

LEGACY_RUN = ref_path("mr2019_input").parents[1]
TAG = "SOUTHAMERICA_GOSAT-BRAZIL_column_20130101"
LEVELS = 17
EDGES = "nesw"


def read_run_fields(run):
    """The fields file of one NAME 8.7.1 run (exactly one is expected)."""
    files = sorted((run / "output").glob("Fields_grid1_C1_T1_*.txt"))
    if len(files) != 1:
        raise ValueError(f"Expected one fields file in {run / 'output'}, found {files}")
    return read_fields(files[0])


def process_run(run):
    """Legacy-process one stage 1 run; returns (processed Dataset as written, raw fields)."""
    fields = read_run_fields(run)
    particles = legacy_particles_from_fate(
        *read_name_fate_file(run / "output" / "Particles_C1.txt")
    )
    met = read_legacy_met_folder(LEGACY_RUN / "Met" / f"Met_{TAG}")
    obs = xr.open_dataset(ref_path("obs")).load()
    result = process_legacy_inputs(fields, particles, met, obs, release_seconds=60.0)
    path = run / "fp_20130101_legacy_processing.nc"
    write_legacy_fp_file(result, path)
    # Compare the file as written: it has the same encoding (rounding) as the production file.
    return xr.open_dataset(path).load(), fields


def arrays(ds):
    """Footprints (sounding, lat, lon) and edge histograms {edge: (sounding, position, height)}."""
    fp = ds["fp"].transpose("time", "lat", "lon").values
    hist = {}
    for e in EDGES:
        da = ds[f"particle_locations_{e}"]
        along = "lon" if e in "ns" else "lat"
        hist[e] = da.transpose("time", along, "height").values
    return fp, hist


def per_sounding_table(A, B, production, legacy_sources):
    fa, ha = arrays(A)
    fb, hb = arrays(B)
    fp, hp = arrays(production)
    lat, lon = B["lat"].values, B["lon"].values
    s_lat, s_lon = B["release_lat"].values, B["release_lon"].values
    binned = ~np.isclose(legacy_sources["dy"].to_numpy()[::LEVELS], 0.0944287686215)
    release_end = pd.to_datetime(legacy_sources["stop_time"].to_numpy()[::LEVELS])
    return pd.DataFrame(
        {
            "release_end": release_end.strftime("%H:%M:%S"),
            "binned": binned,
            "A/B total": total_ratio(fa, fb),
            "A/B near": near_field_ratio(fa, fb, lat, lon, s_lat, s_lon),
            "B/prod total": total_ratio(fb, fp),
            "B/prod near": near_field_ratio(fb, fp, lat, lon, s_lat, s_lon),
            "B-prod shape": shape_difference(fb, fp),
            "hist A-B": histogram_difference(ha, hb),
            "hist B-prod": histogram_difference(hb, hp),
        },
        index=pd.RangeIndex(1, B.sizes["time"] + 1, name="sounding"),
    )


def noise_test(B, B2, production, fields_b, fields_b2, legacy_sources):
    """Print B - B2 (pure noise) beside B - production for every measure."""
    fb, hb = arrays(B)
    f2, h2 = arrays(B2)
    fp, hp = arrays(production)

    print("\n=== Noise test: B vs B2 (seed only) beside B vs production ===")
    print("Per-sounding column totals, std of log ratio:")
    print(f"  B/B2 {np.log(total_ratio(fb, f2)).std():.3f}   "
          f"B/prod {np.log(total_ratio(fb, fp)).std():.3f}")  # fmt: skip

    print("Footprint shape difference, median over soundings (0 = identical, 2 = no overlap):")
    for factor in [1, 2, 4, 8, 16]:
        noise = np.median(shape_difference(fb, f2, factor))
        versus = np.median(shape_difference(fb, fp, factor))
        print(f"  cells {factor:2d}x bigger: B-B2 {noise:.3f}   B-prod {versus:.3f}")

    print("Edge histograms, median difference:")
    print(f"  B-B2 {np.median(histogram_difference(hb, h2)):.3f}   "
          f"B-prod {np.median(histogram_difference(hb, hp)):.3f}")  # fmt: skip

    # Per release, from the raw fields (sources in the same order in all three runs).
    production_fields = read_fields(ref_path("mr2019_fields"))
    totals = pd.DataFrame(
        {
            "level": np.tile(np.arange(1, LEVELS + 1), len(legacy_sources) // LEVELS),
            "rate": legacy_sources["particles_per_second"].to_numpy(),
            "B": fields_b.sum(("lat", "lon")).values,
            "B2": fields_b2.sum(("lat", "lon")).values,
            "prod": production_fields.sum(("lat", "lon")).values,
        }
    )
    rows = []
    for level, group in totals.groupby("level"):
        both = group[(group["B"] > 0) & (group["B2"] > 0) & (group["prod"] > 0)]
        noise = np.log(both["B"] / both["B2"])
        versus = np.log(both["B"] / both["prod"])
        rows.append(
            {
                "level": level,
                "n": len(both),
                "std log B/B2": noise.std(),
                "std log B/prod": versus.std(),
                "mean log B/prod": versus.mean(),
                "sem": versus.std() / np.sqrt(len(both)),
                "zero B/B2/prod": f"{(group['B'] == 0).sum()}/{(group['B2'] == 0).sum()}/"
                f"{(group['prod'] == 0).sum()}",
            }
        )
    by_level = pd.DataFrame(rows).set_index("level")
    print("Per-release footprint totals by level:")
    print(by_level.round(3).to_string())
    for rate, group in totals.groupby("rate"):
        print(f"  {rate}/s releases: ΣB/Σprod {group['B'].sum() / group['prod'].sum():.3f}, "
              f"ΣB2/Σprod {group['B2'].sum() / group['prod'].sum():.3f}, "
              f"ΣB/ΣB2 {group['B'].sum() / group['B2'].sum():.3f}")  # fmt: skip
    return by_level


def main(work):
    work = Path(work)
    A, _ = process_run(work / "run_A")
    B, fields_b = process_run(work / "run_B")
    production = xr.open_dataset(ref_path("production_201301")).sel(time=B["time"].values)
    legacy_sources = read_name_input_sources(ref_path("mr2019_input"))

    table = per_sounding_table(A, B, production, legacy_sources)
    pd.set_option("display.width", 200)
    print(table.round(3).to_string())
    table.to_csv(work / "stage1_comparison.csv")

    for name in ["pressure", "temperature", "PBLH"]:
        same = np.array_equal(B[name].values, production[name].values)
        print(f"{name}: B == production: {same}")

    if not list((work / "run_B2" / "output").glob("Fields_grid1_*.txt")):
        print("\nrun_B2 has no output yet; noise test skipped")
        return
    B2, fields_b2 = process_run(work / "run_B2")
    by_level = noise_test(B, B2, production, fields_b, fields_b2, legacy_sources)
    by_level.to_csv(work / "stage1_noise_by_level.csv")


if __name__ == "__main__":
    main(sys.argv[1])
