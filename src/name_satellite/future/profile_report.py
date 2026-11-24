"""Profile of future-method NAME runs: where the time goes (compute vs met reading vs output).

Where this sits: after the runs. Reads each run folder's ``profile.csv`` (``future.sampler``) and
``inputTiming.txt`` (NAME's timer summary), and writes ``profile_summary.csv`` and one timeline
PNG per run into OUT_FOLDER::

    uv run --with matplotlib python -m name_satellite.future.profile_report MONTH_FOLDER OUT_FOLDER

How the time is split (from the samples, every 5 s): ``busy_threads`` = CPU-seconds per second.

* compute: at least 75% of the run's threads busy (particles moving in parallel);
* serial reading: fewer than 1.5 threads busy while NAME reads data (``rchar`` rising): one
  thread reads met while the others wait;
* serial other: fewer than 1.5 threads busy and no reading (e.g. waiting on the file system, or
  serial work);
* partial: in between.

NAME's own timers (``Use Timers? = Yes``) give the particle loop, output processing and the total;
with the default on-demand met method its met timers stay at zero, which is why the samples are
needed. ``untimed`` = total - particle loop - output processing.
"""

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TIMERS = {
    "TotalTime": "name_total_s",
    "WorkerThread_LPPATs": "particle_loop_s",
    "WorkerThread_ProcessAndOutputResults": "output_s",
}


def read_timers(path):
    """{timer name: total seconds} from NAME's timer summary."""
    values = {}
    current = None
    for line in Path(path).read_text().splitlines():
        match = re.match(r"\s*Timer (\S+) summary:", line)
        if match:
            current = match.group(1)
        elif current and line.strip().startswith("Total time"):
            values[current] = float(line.split()[-1])
            current = None
    return values


def summarise_run(folder):
    """One row of numbers for a run folder with profile.csv (see the module docstring)."""
    folder = Path(folder)
    info = json.loads((folder / "day.json").read_text())
    threads = info["options"]["threads"]
    p = pd.read_csv(folder / "profile.csv")
    dt = p["wall_s"].diff().fillna(0.0)
    reading = p["rchar"].diff().fillna(0) > 1e6
    busy = p["busy_threads"]
    compute = busy >= 0.75 * threads
    serial = busy < 1.5
    row = {
        "run": folder.name,
        "threads": threads,
        "wall_min": p["wall_s"].iloc[-1] / 60,
        "cpu_min": p["cpu_s"].iloc[-1] / 60,
        "mean_busy_threads": p["cpu_s"].iloc[-1] / max(p["wall_s"].iloc[-1], 1e-6),
        "share_compute": dt[compute].sum() / dt.sum(),
        "share_serial_reading": dt[serial & reading].sum() / dt.sum(),
        "share_serial_other": dt[serial & ~reading].sum() / dt.sum(),
        "share_partial": dt[~compute & ~serial].sum() / dt.sum(),
        "read_gb": p["rchar"].iloc[-1] / 1e9,
        "written_gb": p["wchar"].iloc[-1] / 1e9,
        "read_gb_per_min_while_reading": (p["rchar"].diff()[reading].sum() / 1e9)
        / max(dt[reading].sum() / 60, 1e-6),
        "met_files_read": int(p["met_reads"].iloc[-1]),
        "peak_memory_gb": p["hwm_kb"].max() / 1e6,
        "changes": json.dumps(info.get("changes", {})),
    }
    timing = folder / "inputTiming.txt"
    if timing.exists():
        timers = read_timers(timing)
        for name, column in TIMERS.items():
            row[column] = timers.get(name, np.nan)
        row["untimed_s"] = row["name_total_s"] - row["particle_loop_s"] - row["output_s"]
    return row


def timeline(folder, out):
    """Busy threads, GB read and met files read against time for one run."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p = pd.read_csv(Path(folder) / "profile.csv")
    minutes = p["wall_s"] / 60
    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    axes[0].plot(minutes, p["busy_threads"], color="#2a78d6", linewidth=1)
    axes[0].set_ylabel("Busy threads")
    axes[1].plot(minutes, p["rchar"] / 1e9, color="#eb6834")
    axes[1].set_ylabel("Read (GB)")
    axes[2].plot(minutes, p["met_reads"], color="#1f9e74")
    axes[2].set_ylabel("Met files read")
    axes[2].set_xlabel("Minutes since NAME started")
    axes[0].set_title(Path(folder).name, loc="left")
    for ax in axes:
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(Path(out) / f"timeline_{Path(folder).name}.png", dpi=100)
    plt.close(fig)


def main(month_folder, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    folders = sorted(f for f in Path(month_folder).iterdir() if (f / "profile.csv").exists())
    rows = []
    for folder in folders:
        if len(pd.read_csv(folder / "profile.csv")) < 3:
            print(f"{folder.name}: fewer than 3 samples (run failed or too short); skipped")
            continue
        rows.append(summarise_run(folder))
        timeline(folder, out)
    table = pd.DataFrame(rows)
    table.to_csv(out / "profile_summary.csv", index=False)
    print(table.round(3).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
