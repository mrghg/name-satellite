"""Write two short NAME test runs for the GLIDE reference runs (debug queue).

1. ``identical``: two identical point sources (same place, height and hour) at Mace Head in one
   ``Fixed`` run, 6 h of travel. GLIDE spec §10 check 6: if their footprints are identical, the
   sources share random numbers in a ``Fixed`` run, and `main` must use ``Random`` (§9).
2. ``column``: one source at Lamont spanning 950-100 hPa (Z-coord Pa), released for 5 minutes at
   1,000/s, with every particle's position (in Pa) at the end of the release. The histogram of
   starting pressures shows whether NAME spreads particles evenly in pressure or in height
   (§5: one source, or 28 stacked layers). 950 hPa is safely above Lamont's model surface
   (about 970 hPa); this test is about the spread, not the exact surface pressure.

    uv run python dev/scripts/glide_tests/prepare_source_tests.py identical RUN_FOLDER
    uv run python dev/scripts/glide_tests/prepare_source_tests.py column RUN_FOLDER
    uv run python dev/scripts/glide_tests/prepare_source_tests.py memory RUN_FOLDER

3. ``memory``: 24 hourly releases at Mace Head at 10x the spec's rate (60/s), 24 h of travel, on
   8 threads (QoS high). Near its end about 24 x 216,000 = 5.2 M particles are alive, as many as a
   block run of one point receptor (240 live releases x 21,600). Measures memory (MaxRSS) and
   run time with the spec's settings (deep convection New), to size `Max # Particles` and jobs.

Then in the folder: ``r=$(sbatch --parsable restore_met.sbatch)``, ``sbatch
--dependency=afterok:$r job.sbatch``. Results and how they were read:
dev/notes/name-behaviour.md and docs/glide-reference-runs.md (16 Nov 2026).
"""

import sys
from dataclasses import replace

import pandas as pd

from name_satellite.glide_runs import spec
from name_satellite.glide_runs.plan import point_release_plan
from name_satellite.glide_runs.prepare import write_run_folder
from name_satellite.name_input import ParticlePositionOutput


def identical_sources(run_folder):
    receptor = spec.RECEPTORS["MHD"]
    one = point_release_plan(receptor, [pd.Timestamp("2019-01-08 12:00")])
    two = one.copy()
    two["source_name"] = "SourceID_2"
    plan = pd.concat([one, two], ignore_index=True)
    settings = spec.reference_settings(
        run_name="identical",
        output_folder=f"{run_folder}/output",
        domain=spec.EUR,
        met_folder=f"{run_folder}/met",
        max_particles=spec.MAX_PARTICLES_POINT,
    )
    settings = replace(settings, max_age_hours=6)
    write_run_folder(
        run_folder, plan, settings, "glide_identical", "debug", "00:50:00", "16G",
        record={"test": "two identical sources, Fixed seed (spec §10 check 6)"},
    )  # fmt: skip


def column_spread(run_folder):
    receptor = spec.RECEPTORS["SGP"]
    start = pd.Timestamp("2019-07-08 19:00")
    bottom_pa, top_pa = 95000.0, 10000.0
    plan = pd.DataFrame(
        {
            "source_name": ["SourceID_1"],
            "release_time": [start],
            "release_start": [start + pd.Timedelta(minutes=5)],  # back-run: the later time
            "release_end": [start],
            "receptor": [receptor.code],
            "lon": [receptor.lon],
            "lat": [receptor.lat],
            "z_coord": ["Pa"],
            "z": [(bottom_pa + top_pa) / 2],
            "dz": [bottom_pa - top_pa],
            "dlon": [0.0],
            "dlat": [0.0],
            "particles_per_second": [1000],
        }
    )
    settings = spec.reference_settings(
        run_name="column",
        output_folder=f"{run_folder}/output",
        domain=spec.NAM,
        met_folder=f"{run_folder}/met",
        max_particles=spec.MAX_PARTICLES_POINT,
    )
    # Positions in Pa at the end of the release (19:00) and one sync step later (18:55).
    positions = ParticlePositionOutput(start, n_times=2, interval_minutes=5, z_coord="Pa")
    settings = replace(settings, max_age_hours=1, particle_positions=positions, fate_met=True)
    write_run_folder(
        run_folder, plan, settings, "glide_column", "debug", "00:50:00", "16G",
        record={"test": "spread of a 950-100 hPa source (spec §5)"},
    )  # fmt: skip


def memory_probe(run_folder):
    times = pd.date_range("2019-01-08 00:00", "2019-01-08 23:00", freq="1h")
    plan = point_release_plan(spec.RECEPTORS["MHD"], times)
    plan["particles_per_second"] = 10 * spec.POINT_PARTICLES_PER_SECOND
    settings = spec.reference_settings(
        run_name="memory",
        output_folder=f"{run_folder}/output",
        domain=spec.EUR,
        met_folder=f"{run_folder}/met",
        max_particles=spec.MAX_PARTICLES_POINT,
        threads=8,
    )
    settings = replace(settings, max_age_hours=24)
    write_run_folder(
        run_folder, plan, settings, "glide_memory", "high", "06:00:00", "96G",
        record={"test": "memory probe: ~5.2 M live particles (24 releases at 60/s, 24 h)"},
    )  # fmt: skip


if __name__ == "__main__":
    test, folder = sys.argv[1:3]
    tests = {"identical": identical_sources, "column": column_spread, "memory": memory_probe}
    tests[test](folder)
    print(f"Wrote {folder}")
