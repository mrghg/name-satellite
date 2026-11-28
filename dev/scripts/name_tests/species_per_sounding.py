"""Species per sounding: can NAME weight the levels itself, with far fewer field requests?

(NAME test, future-method branch.) A future-method run has one species and, per release, three
field requests restricted to that release (its group's 0.1° and 0.5° boxes, the global 2.5°
grid); each sounding's column footprint is formed afterwards as Σ over levels of pressure weight
x averaging kernel x the level's footprint. Run time grows with particles x field requests in the
run. Here instead each sounding has its own species: each release of a group emits the species of
every sounding sharing that group, at a rate equal to that sounding's weight for the release's
level (g/s), and there are three requests per **sounding** (species, no source restriction). The
field of a sounding's species should then be the weighted sum directly.

Writes WORK/<layout>/{input.txt, job.sbatch} for the layouts:

* ``sources``: the future method as it is (requests per release, one species);
* ``species``: one species per sounding, requests per sounding.

Both use the same releases, particles and random seed, so the particles' paths are the same and
the two should give the same column fields (to the precision of NAME's text output).

    uv run python dev/scripts/name_tests/species_per_sounding.py PART_FOLDER WORK [--hours=1]

PART_FOLDER: a completed future-method run folder (its ``release_plan.csv``, ``grids.csv`` and
``day.json``), e.g. a ``profile_split4_<day>_part<k>``. ``--hours``: travel time (1: debug queue).
"""

import dataclasses
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from name_satellite.domains import SITES
from name_satellite.future.setup import (
    PATCHED_NAME_DIR,
    Options,
    input_text,
    run_settings,
)
from name_satellite.lotus import write_job_script
from name_satellite.releases import MODELLED_LEVELS


def sounding_weights(part_folder, info, plan):
    """Weights (pressure weight x averaging kernel) by sounding and level, for the soundings whose
    release groups are in the plan. Returns a DataFrame (sounding_id, release_group, level, weight).
    """
    day_folder = Path(info["day_folders"][0])
    day = pd.Timestamp(info["day"])
    members = pd.read_csv(day_folder / "release_groups.csv", dtype=str)
    members = members[members["release_group"].isin(set(plan["release_group"]))]
    with xr.open_dataset(day_folder / f"obs_{info['site']}_{day:%Y%m%d}.nc") as f:
        obs = f.load()
    n_levels = MODELLED_LEVELS[info["product"]]
    ids = [str(s) for s in obs["sounding_id"].values]
    weights = (obs["pressure_weights"].values[:, :n_levels]
               * obs["xch4_averaging_kernel"].values[:, :n_levels]).astype(float)  # fmt: skip
    rows = []
    for sounding, group in zip(members["sounding_id"], members["release_group"], strict=True):
        for level in range(1, n_levels + 1):
            w = weights[ids.index(sounding), level - 1]
            if not np.isfinite(w) or w < 0:
                raise ValueError(f"Weight {w} for sounding {sounding} level {level}")
            rows.append((sounding, group, level, w))
    return pd.DataFrame(rows, columns=["sounding_id", "release_group", "level_number", "weight"])


def species_layout(text, plan, boxes, weights):
    """Rewrite a future-method input (one species, requests per release) as species per sounding."""
    lines = text.split("\n")
    soundings = list(dict.fromkeys(weights["sounding_id"]))
    name_of = {s: f"SND{k + 1}" for k, s in enumerate(soundings)}

    # Main Options: room for every species.
    i = lines.index(next(line for line in lines if line.startswith("Absolute or Relative Time?")))
    lines[i] = lines[i].replace("Max # Sources,", "Max # Species, Max # Sources,")
    head = lines[i].split(", ")
    values = lines[i + 1].split(", ")
    values.insert(head.index("Max # Species"), str(len(soundings)))
    lines[i + 1] = ", ".join(values)

    # Species and Species Uses: one inert tracer per sounding (as INERT-TRACER).
    i = lines.index("INERT-TRACER, INERT, Stable, 0.00E+00, , 0.0, 1.0, g")
    lines[i : i + 1] = [f"{name_of[s]}, INERT, Stable, 0.00E+00, , 0.0, 1.0, g" for s in soundings]
    i = lines.index("INERT-TRACER, Yes, No, No")
    lines[i : i + 1] = [f"{name_of[s]}, Yes, No, No" for s in soundings]

    # Sources: one Source Strength column per sounding sharing the release's group.
    per_group = weights.groupby("release_group")["sounding_id"].unique()
    n_columns = max(len(v) for v in per_group)
    i = next(k for k, line in enumerate(lines) if line.startswith("Name, Start Time, Stop Time"))
    lines[i] = lines[i].replace(
        "Source Strength,", ", ".join(["Source Strength"] * n_columns) + ","
    )
    w = weights.set_index(["sounding_id", "level_number"])["weight"]
    for release in plan.itertuples():
        k = next(k for k, line in enumerate(lines) if line.startswith(f"{release.source_name}, "))
        strengths = [f"{name_of[s]} {w.loc[(s, release.level_number)]:.6e} g/s"
                     for s in per_group[release.release_group]]  # fmt: skip
        strengths += [""] * (n_columns - len(strengths))
        if lines[k].count("INERT-TRACER 1.0 g/s") != 1:
            raise ValueError(f"Unexpected source line for {release.source_name}")
        lines[k] = lines[k].replace("INERT-TRACER 1.0 g/s", ", ".join(strengths))

    # Field requests: three per sounding (its group's boxes and the global grid), by species.
    requests = [k for k, line in enumerate(lines) if line.startswith("Req_")]
    template = lines[requests[0]].split(", ")
    new = []
    by_group = boxes.set_index("release_group")
    for s in soundings:
        group = by_group.loc[weights.loc[weights["sounding_id"] == s, "release_group"].iloc[0]]
        for grid in (group["fine_grid"], group["medium_grid"], "HGridC"):
            fields = list(template)
            fields[0], fields[2], fields[3], fields[4] = (
                f"Req_{grid}_{name_of[s]}",
                name_of[s],
                "",
                grid,
            )
            fields[-1] = f"Fields_{grid}"
            new.append(", ".join(fields))
    lines[requests[0] : requests[-1] + 1] = new
    if requests != list(range(requests[0], requests[-1] + 1)):
        raise ValueError("Field requests are not one consecutive block")
    return "\n".join(lines), name_of


def main(part_folder, work, hours):
    part_folder, work = Path(part_folder), Path(work)
    info = json.loads((part_folder / "day.json").read_text())
    plan = pd.read_csv(part_folder / "release_plan.csv", dtype={"release_group": str},
                       parse_dates=["release_start", "release_end"])  # fmt: skip
    boxes = pd.read_csv(part_folder / "grids.csv", dtype={"release_group": str})
    weights = sounding_weights(part_folder, info, plan)
    # Both layouts start from the "sources" input; species_layout rewrites it. (future.setup now
    # writes the species layout itself: Options.layout, the default since 27 Nov 2026.)
    options = dataclasses.replace(Options(**info["options"]), test=hours == 1, layout="sources")
    site = SITES[info["site"]]
    day = pd.Timestamp(info["day"])
    for layout in ("sources", "species"):
        run_folder = work / layout
        (run_folder / "output").mkdir(parents=True, exist_ok=False)
        settings = run_settings(site, run_folder, day, day, options)
        if hours != 1:
            settings = dataclasses.replace(settings, max_age_hours=hours)
        text = input_text(plan, boxes, settings, options)
        if layout == "species":
            text, name_of = species_layout(text, plan, boxes, weights)
            pd.Series(name_of, name="species").rename_axis("sounding_id").to_csv(
                run_folder / "species.csv"
            )
        (run_folder / "input.txt").write_text(text)
        n_requests = len(re.findall(r"^Req_", text, flags=re.M))
        qos, time = ("debug", "00:30:00") if hours == 1 else ("high", "24:00:00")
        write_job_script(run_folder / "job.sbatch", run_folder, run_folder / "input.txt",
                         threads=options.threads, qos=qos, time=time, memory="48G",
                         job_name=f"spec_{layout}", name_dir=PATCHED_NAME_DIR)  # fmt: skip
        print(f"{run_folder}: {len(plan)} releases, {n_requests} field requests")
    weights.to_csv(work / "weights.csv", index=False)


if __name__ == "__main__":
    flags = dict(a.removeprefix("--").split("=", 1) for a in sys.argv[3:] if "=" in a)
    main(sys.argv[1], sys.argv[2], int(flags.get("hours", 1)))
