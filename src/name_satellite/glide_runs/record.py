"""The NAME run record for GLIDE (spec §8, "Run record"): one YAML file.

Where this sits
---------------
    prepare.py -> LOTUS2 -> deliver.py  ->  **record.py: what was run, how, and what NAME saw**

The record contains NAME's complete input file, which includes the NAME team's met declarations
(licensed with NAME). It is therefore written next to the NAME output on the ACRG group workspace
and never committed to a git repository (spec §8).

What goes in (spec §8 list): model and patch; one complete input file; the met (series, grid,
levels, interval); the receptor's model ground height and surface type; the Lamont surface
pressures (if any); the domain top and what happens there; the main settings; the unit
conversion; the random seed setting of each run; wall time, threads and hardware of each run
(from Slurm, ``sacct``); departures from the spec. Values measured from the met files themselves
are passed in (``MET_DESCRIPTION``) with the date they were read.
"""

import json
import subprocess
from pathlib import Path

import pandas as pd
import yaml

from name_satellite.glide_runs import spec
from name_satellite.glide_runs.column import COLUMN_METHOD
from name_satellite.release_met import read_release_met_file

# Read from the restored Mk10 files and NAME 8.7.1's met definition on 16 Nov 2026 (PP header of
# MO201812281800.UMG_Mk10_I_L59PT1.pp; file times every 3 h).
MET_DESCRIPTION = {
    "source": "Met Office Unified Model global analyses, NAME met archive on JASMIN",
    "series": "UMG_Mk10PT (one series for every block)",
    "grid": "regular latitude-longitude, 2,560 points around the globe; 0.140625 deg (lon) x "
    "0.09375 deg (lat)",
    "levels": 59,
    "interval_between_fields_hours": 3,
}

PATCH_DESCRIPTION = (
    "ACRG patch to NAME 8.7.1: adds each particle's per-source index and its stop reason (left "
    "the domain / reached the maximum travel time) to the particle-fate output. It does not "
    "change dispersion: on the one-day check (16 Nov 2026) the patched and stock builds wrote "
    "identical fate records for all 518,424 particles, and identical footprints except one value "
    "differing in the 7th significant figure (10.50166 vs 10.50167 ppm s)."
)

TOP_DESCRIPTION = (
    "NAME's computational domain top is 20,000 m above the model's ground. A particle that "
    "rises above it is removed; its end record has fate 'top'."
)

CONVERSION = (
    "srr = NAME time-integrated mixing ratio in the 0-40 m layer (ppm s; source 1 g/s of a "
    "tracer of molar mass 1 g/mol, i.e. 1 mol/s) x 1e-6 x cell area (m2, sphere of radius "
    "6,367,500 m) / 3,600 s. NAME divides each source's mass among the particles it actually "
    "released, so no further normalisation by particle number is applied."
)


def slurm_usage(job_ids):
    """Wall time, CPUs, node and peak memory of Slurm jobs (``sacct``), as a list of dicts."""
    out = subprocess.run(
        [
            "sacct", "-n", "-P", "-j", ",".join(str(j) for j in job_ids),
            "-o", "JobID,JobName,Elapsed,AllocCPUS,NodeList,MaxRSS,State",
        ],
        capture_output=True, text=True, check=True,
    ).stdout  # fmt: skip
    rows = {}
    peak = {}
    for line in out.splitlines():
        job, name, elapsed, cpus, node, rss, state = line.split("|")
        if "." in job:  # a job step (.batch, .extern, .0): NAME's memory is in its step
            if rss:
                peak[job.split(".")[0]] = max(peak.get(job.split(".")[0], ""), rss, key=_kilobytes)
            continue
        rows[job] = {"job": job, "name": name, "elapsed": elapsed, "cpus": cpus, "node": node,
                     "state": state}  # fmt: skip
    for job, row in rows.items():
        row["max_rss"] = peak.get(job, "")
    return list(rows.values())


def _kilobytes(rss):
    """Slurm's MaxRSS text (e.g. "5023852K", "18M") in kB, for comparing."""
    if not rss:
        return 0.0
    scale = {"K": 1.0, "M": 1024.0, "G": 1024.0**2}.get(rss[-1], 1.0 / 1024.0)
    return float(rss.rstrip("KMG")) * scale


def mesoscale_parameters(definition_file):
    """NAME's unresolved-mesoscale parameters for a met series, read from its met definition.

    Read at run time from the NAME installation (never copied into this repository): the
    "Mesoscale SigU" (m/s) and "Mesoscale TauU" (s) columns of the definition's met-file rows.
    Returns a short text; fails if the rows disagree.
    """
    path = Path(spec.MET_VERSION_DEFINITIONS) / definition_file
    lines = path.read_text().splitlines()
    heading = next(i for i, line in enumerate(lines) if "Mesoscale SigU" in line)
    names = [c.strip() for c in lines[heading].split(",")]
    values = set()
    for line in lines[heading + 1 :]:
        if not line.strip():
            break
        cells = [c.strip() for c in line.split(",")]
        values.add((cells[names.index("Mesoscale SigU")], cells[names.index("Mesoscale TauU")]))
    if len(values) != 1:
        raise ValueError(f"Mesoscale parameters differ between met files in {path.name}")
    sig_u, tau_u = values.pop()
    return f"velocity standard deviation {sig_u} m/s, Lagrangian timescale {tau_u} s"


def receptor_met(run_folder):
    """Model ground height (m above sea level) and land fraction at the receptor, from the run."""
    columns, values = read_release_met_file(Path(run_folder) / "output" / "ReleaseMet_R1_C1.txt")
    row = dict(zip(columns["request"], values.loc["R1"], strict=True))
    land = row["R1_LAND"]
    return {
        "model_ground_height_m_asl": float(row["R1_TOPO"]),
        "surface_type": (
            "not available: NAME reports Land Fraction as missing with Mk10 met"
            if land < -1e30
            else f"land fraction {land}"
        ),
    }


def build_record(runs, departures, notes=None, surface_pressure_csv=None):
    """The run record as a dict.

    Args:
        runs: list of dicts, one per run: ``folder`` (run folder), ``job`` (Slurm job id) and
            ``purpose`` (e.g. "main", "§10 check 5: second Fixed run").
        departures: list of strings, each a departure from the spec with its reason.
        notes: other entries (dict), e.g. the §10 check results.
        surface_pressure_csv: Lamont surface pressures (``column.read_surface_pressure_runs``).
    """
    first = Path(runs[0]["folder"])
    first_record = json.loads((first / "run.json").read_text())
    settings = first_record["settings"]
    meso = mesoscale_parameters(spec.MET_VERSION.definition_file)
    receptors = {}
    run_entries = []
    usage = {row["job"]: row for row in slurm_usage([r["job"] for r in runs])}
    for run in runs:
        folder = Path(run["folder"])
        record = json.loads((folder / "run.json").read_text())
        code = record["receptor"]["code"]
        receptors.setdefault(code, receptor_met(folder))
        run_entries.append(
            {
                "folder": str(folder),
                "purpose": run["purpose"],
                "receptor": code,
                "releases": f"{record['first_release']} to {record['last_release']} "
                f"({record['n_releases']})",
                "variant": record["variant"],
                "random_seed": record["settings"]["random_seed"],
                "name_build": record["name_build"],
                "slurm": usage.get(str(run["job"]), {}),
                "hardware": "JASMIN LOTUS2 (node in slurm.node)",
            }
        )
    result = {
        "specification": spec.SPEC_VERSION,
        "model": {
            "name": "NAME III",
            "version": "8.7.1",
            "patch": PATCH_DESCRIPTION,
            "executable": "nameiii_64bit_par.exe (OpenMP)",
        },
        "meteorology": MET_DESCRIPTION,
        "receptors": receptors,
        "domain_top": TOP_DESCRIPTION,
        "settings": {
            "sync_time_s": int(settings["sync_minutes"]) * 60,
            "max_travel_time_h": int(settings["max_age_hours"]),
            "boundary_layer_depth_limits_m": [
                settings["min_bl_depth_m"],
                settings["max_bl_depth_m"],
            ],
            "boundary_layer_depth_source": "UM field (Use NWP BL Depth? = Yes), then limited",
            "deep_convection": f"{settings['deep_convection']} (NAME's mass-flux scheme)",
            "turbulence": "on (NAME's standard scheme)",
            "unresolved_mesoscale_motions": (
                f"on; NAME's standard parameters for this met (its met definition): {meso}"
                if settings["mesoscale_motions"]
                else "off"
            ),
            "max_particles": settings["max_particles"],
            "particles_per_second": spec.POINT_PARTICLES_PER_SECOND,
            "particles_per_second_column": spec.COLUMN_PARTICLES_PER_SECOND,
            "output_layer": "0-40 m above the model's ground, one layer",
            "output_grid": "spec §2 (EUR 279 x 159, NAM 259 x 143, 0.25 deg)",
            "column_method": COLUMN_METHOD,
        },
        "unit_conversion": CONVERSION,
        "runs": run_entries,
        "departures_from_spec": departures,
        "input_file_of_first_run": (first / "input.txt").read_text(),
    }
    if surface_pressure_csv is not None:
        table = pd.read_csv(surface_pressure_csv)
        result["lamont_surface_pressure_hpa"] = {
            row.release_time: round(row.surface_pressure_pa / 100.0, 3)
            for row in table.itertuples()
        }
    if notes:
        result["notes"] = notes
    return result


def write_record(path, record):
    """Write ``record`` as YAML (long strings as literal blocks, for the input file)."""

    def literal(dumper, text):
        style = "|" if "\n" in text else None
        return dumper.represent_scalar("tag:yaml.org,2002:str", text, style=style)

    yaml.add_representer(str, literal)
    Path(path).write_text(yaml.dump(record, sort_keys=False, allow_unicode=True, width=100))
    return Path(path)
