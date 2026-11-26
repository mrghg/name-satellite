"""How many field requests can a NAME run hold? (NAME test, future-method branch)

A 3-day batch (1,520 releases, three requests each: its group's 0.1° and 0.5° boxes and the global
2.5° grid; 4,560 in all) stopped with a fatal error saying a field-request limit "is not large
enough"; single days with up to 944 releases (2,832 requests) ran. This writes 1-hour test runs
(debug queue) built from that batch's releases, one per CASE:

* ``<N>``: the first N releases' coarse requests on HGridC1, the rest on an identical copy HGridC2
  (own output group). Does the limit count requests per horizontal grid?
* ``<N>split``: all coarse requests on HGridC1, the first N in one output group and the rest in a
  second. Per output group?
* ``total<T>``: exactly T requests in the run (whole releases' three requests, then coarse-only
  requests to make up T; other releases get none). Per run?
* ``tgrid``: all 4,560 requests, the coarse ones on a copy of the time grid (TGridEndOfRun2).
  Per time grid?

    uv run python dev/scripts/name_tests/field_request_limit.py BATCH_FOLDER WORK CASE [CASE ...]

BATCH_FOLDER: the month's ``profile_batch3days`` folder (its day folders hold the surface-pressure
runs). Writes WORK/N<CASE>/{input.txt, job.sbatch}; submit with ``sbatch``.

Results (26 Nov 2026, patched 8.7.1): see dev/notes/name-behaviour.md, "Limits met in our runs".
"""

import json
import re
import sys
from pathlib import Path

import pandas as pd

from name_satellite.batching import combine_release_plans
from name_satellite.domains import SITES
from name_satellite.future.grids import COARSE, global_box, name_grid_line
from name_satellite.future.setup import (
    PATCHED_NAME_DIR,
    Options,
    final_plan,
    group_boxes,
    input_text,
    run_settings,
)
from name_satellite.lotus import write_job_script


def batch_input_lines(batch_folder, run_folder):
    """The batch's NAME input with 1 h of travel and every coarse request on HGridC1, as lines;
    and the number of releases."""
    info = json.loads((batch_folder / "day.json").read_text())
    site = SITES[info["site"]]
    options = Options(**{**info["options"], "test": True})
    day_folders = [Path(p) for p in info["day_folders"]]
    plans = {f.name: final_plan(f, site, info["product"], write_tables=False) for f in day_folders}
    plan = combine_release_plans(plans)
    boxes = group_boxes(plan)
    boxes["coarse_grid"] = "HGridC1"
    days = [pd.Timestamp(json.loads((f / "day.json").read_text())["day"]) for f in day_folders]
    settings = run_settings(site, run_folder, min(days), max(days), options)
    return input_text(plan, boxes, settings, options).split("\n"), len(plan)


def request_lines(lines, prefix):
    """Line numbers of the field requests whose name starts with ``prefix``."""
    return [i for i, line in enumerate(lines) if line.startswith(prefix)]


def build_input(batch_folder, run_folder, case):
    """Input text for one CASE (see the module docstring) and a one-line description."""
    lines, n_releases = batch_input_lines(batch_folder, run_folder)
    coarse = request_lines(lines, "Req_HGridC1_")
    if len(coarse) != n_releases:
        raise ValueError("Expected one coarse request per release")

    if case == "tgrid":
        for i in coarse:
            lines[i] = lines[i].replace("TGridEndOfRun,", "TGridEndOfRun2,")
        t_line = next(i for i, line in enumerate(lines) if line.startswith("TGridEndOfRun,"))
        lines.insert(t_line + 1, lines[t_line].replace("TGridEndOfRun,", "TGridEndOfRun2,"))
        about = f"{3 * n_releases} requests; the coarse ones on a second time grid"

    elif case.startswith("total"):
        total = int(case.removeprefix("total"))
        full, extra = divmod(total, 3)
        # Requests are written per release (fine, medium, coarse): keep the first `full` releases'
        # three, then `extra` more releases' coarse request only, and drop every other request.
        fine, medium = request_lines(lines, "Req_HGridF"), request_lines(lines, "Req_HGridM")
        drop = set(fine[full:]) | set(medium[full:]) | set(coarse[full + extra :])
        lines = [line for i, line in enumerate(lines) if i not in drop]
        kept = len(request_lines(lines, "Req_"))
        if kept != total:
            raise ValueError(f"Kept {kept} requests, wanted {total}")
        about = f"{kept} requests in the run"

    elif case.endswith("split"):
        n_first = int(case.removesuffix("split"))
        for i in coarse[n_first:]:
            lines[i] = re.sub(r"Fields_HGridC1\s*$", "Fields_HGridC1b", lines[i])
        about = f"{n_first} of {n_releases} coarse requests in the first of two output groups"

    else:
        n_first = int(case)
        for i in coarse[n_first:]:
            lines[i] = lines[i].replace("HGridC1", "HGridC2")
        grid_line = next(i for i, line in enumerate(lines) if line.startswith("HGridC1,"))
        lines.insert(grid_line + 1, name_grid_line("HGridC2", global_box(), COARSE.cell_deg))
        about = f"{n_first} of {n_releases} coarse requests on HGridC1, the rest on HGridC2"
    return "\n".join(lines), about


def main(batch_folder, work, cases):
    batch_folder, work = Path(batch_folder), Path(work)
    for case in cases:
        run_folder = work / f"N{case}"
        run_folder.mkdir(parents=True, exist_ok=False)
        (run_folder / "output").mkdir()  # NAME does not create it
        text, about = build_input(batch_folder, run_folder, case)
        (run_folder / "input.txt").write_text(text)
        write_job_script(
            run_folder / "job.sbatch", run_folder, run_folder / "input.txt", threads=4,
            qos="debug", time="00:30:00", memory="48G", job_name=f"reqlim_{case}",
            name_dir=PATCHED_NAME_DIR,
        )  # fmt: skip
        print(f"{run_folder}: {about}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
