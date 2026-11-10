"""Restore all met files a period needs into a met folder (run as a LOTUS2 job).

Written into jobs by ``lotus.write_met_restore_job``:

    uv run python -m name_satellite.met_restore UMG_Mk6PT 2012-12-02T14:58 2013-01-01T16:46 FOLDER
"""

import sys

from name_satellite.met import MET_VERSIONS, required_met_files, restore_met_files


def main(version_name, start, end, folder):
    version = next((v for v in MET_VERSIONS if v.name == version_name), None)
    if version is None:
        raise SystemExit(f"Unknown met version {version_name!r}")
    files = required_met_files(version, start, end)
    print(f"{len(files)} met files for {version.name}, {start} to {end}")
    print(restore_met_files(version, files, folder))


if __name__ == "__main__":
    main(*sys.argv[1:5])
