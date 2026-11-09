"""Restore all met files needed for a period into a met folder (run as a LOTUS2 job).

uv run python dev/scripts/restore_met.py UMG_Mk6PT 2012-12-02T14:58 2013-01-01T16:46 <folder>
"""

import sys

from name_satellite.met import MET_VERSIONS, required_met_files, restore_met_files

if __name__ == "__main__":
    version_name, start, end, folder = sys.argv[1:5]
    version = next(v for v in MET_VERSIONS if v.name == version_name)
    files = required_met_files(version, start, end)
    print(f"{len(files)} met files for {version.name}, {start} to {end}")
    print(restore_met_files(version, files, folder))
