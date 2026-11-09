"""Compare release-point met from a release-met test run with the legacy met-run files.

For each sounding, levels 1 and 2 (groups 1 and 2 of prepare_release_met_test.py): differences
in temperature, pressure, wind and boundary-layer depth, and, for releases the legacy met run
found underground (legacy pressure ≠ requested), the requested, legacy and new pressures beside
the new run's surface pressure.

    uv run python dev/scripts/compare_release_met.py WORK_FOLDER YYYYMMDD
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from name_satellite.legacy import read_legacy_met_folder
from name_satellite.refdata import ref_path
from name_satellite.release_met import read_release_met

LEGACY_RUN = ref_path("conc2018_input").parents[1]
RUN_NAME = "SOUTHAMERICA_GOSAT-BRAZIL_column"


def main(work, day):
    work = Path(work)
    soundings = pd.read_csv(work / "soundings.csv", parse_dates=["time"])
    new = read_release_met(work / "output", soundings)
    legacy = read_legacy_met_folder(LEGACY_RUN / "Met" / f"Met_{RUN_NAME}_{day}")
    pd.set_option("display.width", 200)

    for level in (1, 2):
        ours = new[new["group"] == level].reset_index(drop=True)
        theirs = legacy[legacy["level_number"] == level].reset_index(drop=True)
        if not (ours["point_number"] == theirs["point_number"]).all():
            raise ValueError("Soundings don't line up with the legacy met files")
        diff = pd.DataFrame(
            {
                "temperature (K)": ours["temperature_k"] - 273.15 - theirs["temperature_c"],
                "pressure (Pa)": ours["pressure_pa"] - theirs["pressure_pa"],
                "wind speed (m/s)": ours["wind_speed_m_s"] - theirs["wind_speed_m_s"],
                "wind direction (deg)": (
                    (ours["wind_direction_deg"] - theirs["wind_direction_deg"] + 180) % 360 - 180
                ),
                "BL depth (m)": ours["boundary_layer_depth_m"] - theirs["boundary_layer_depth_m"],
            }
        )
        print(f"Level {level}, {len(ours)} soundings: largest |new - legacy|")
        print(diff.abs().max().round(4).to_string())

        underground = np.abs(ours["z_pa"] - theirs["pressure_pa"]) > 1.0
        print(f"  releases legacy found underground: {underground.sum()}")
        if underground.any():
            table = pd.DataFrame(
                {
                    "sounding": ours["point_number"],
                    "requested": ours["z_pa"],
                    "legacy": theirs["pressure_pa"],
                    "new": ours["pressure_pa"],
                    "new surface": ours["surface_pressure_pa"],
                }
            )[underground]
            print(table.round(2).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
