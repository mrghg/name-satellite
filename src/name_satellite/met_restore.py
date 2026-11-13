"""Restore all met files a period needs into a met folder (run as a LOTUS2 job).

Where this sits (see docs/workflow.md)
--------------------------------------
    ``prepare`` writes ``restore_met.sbatch``  ->  **this, on LOTUS2, before any NAME run of the
    month**  ->  the surface-pressure and main NAME runs read the met folder

Written into jobs by ``lotus.write_met_restore_job``::

    uv run python -m name_satellite.met_restore VERSION START END FOLDER
    uv run python -m name_satellite.met_restore UMG_Mk6PT 2012-12-02T14:58 2013-01-01T16:46 FOLDER

VERSION is a ``met.MET_VERSIONS`` name; START and END the first and last time of met needed
(UTC); FOLDER the met folder. It decompresses (or links) every met file NAME will ask for
(``met.required_met_files``, ``met.restore_met_files``) and prints the number of files and the
counts restored / linked / already present (in the job's ``restore-<jobid>.out``). Safe to run
again: files already present with the right size are kept.
"""

import sys

from name_satellite.met import MET_VERSIONS, required_met_files, restore_met_files


def main(version_name, start, end, folder):
    """Restore the met of ``version_name`` from ``start`` to ``end`` into ``folder``.

    Raises:
        SystemExit for an unknown version; errors of ``met.restore_met_files``.
    """
    version = next((v for v in MET_VERSIONS if v.name == version_name), None)
    if version is None:
        raise SystemExit(f"Unknown met version {version_name!r}")
    files = required_met_files(version, start, end)
    print(f"{len(files)} met files for {version.name}, {start} to {end}")
    print(restore_met_files(version, files, folder))


if __name__ == "__main__":
    main(*sys.argv[1:5])
