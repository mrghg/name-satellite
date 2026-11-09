"""Which Met Office global (UM) met data a NAME run uses, and whether it is all in the archive.

Where this sits (see docs/workflow.md)
--------------------------------------
    2. plan releases  ->  **choose met**  ->  3. write NAME inputs

A back-run of one day needs met from (earliest release - maximum particle age) to the latest
release. The UM global model is archived in "versions" (Mk6 ... Mk12); a run must use one version
throughout (the legacy scripts shortened the run near changeovers, docs/deviations.md X2).

The version dates below are those of the NAME team's own 8.7.1 demo scripts
(``~/name/examples/DemoScripts_v8_7_1/AirHistoryBackRuns/BackRuns_OneSite_ByDay.sh``), NOT the
legacy ACRG dates (which stopped at Mk10, X1). File names come from NAME's met definition files,
so the archive check uses exactly the files NAME will ask for.
"""

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

NAME_GWS = Path("/gws/ssde/j25a/name")
NAME_DIR = NAME_GWS / "code" / "NAME-8.7.1"
MET_ARCHIVE = NAME_GWS / "met_archive"
DEMO_UTILITIES = NAME_GWS / "examples" / "DemoScripts_v8_7_1" / "CommonUtilities"


@dataclass(frozen=True)
class MetVersion:
    """One UM global met version, as used by the NAME 8.7.1 demo scripts."""

    name: str  # e.g. "UMG_Mk6PT"
    first_day: str  # first and last day this version is used for (inclusive), YYYY-MM-DD
    last_day: str | None  # None = still current
    definition_file: str  # in NAME's Resources/Defns
    declaration_file: str  # in the demo CommonUtilities/MetDeclarations
    archive_folder: str  # under met_archive


MET_VERSIONS = (
    MetVersion("UMG_Mk6PT", "2010-03-09", "2013-04-29", "MetDefnUMG_Mk6_L59PTpp.txt",
               "Use_UMG_Mk6_L59PTpp.txt", "Global/UMG_Mk6PT"),
    MetVersion("UMG_Mk7PT", "2013-04-30", "2014-07-14", "MetDefnUMG_Mk7_L59PTpp.txt",
               "Use_UMG_Mk7_L59PTpp.txt", "Global/UMG_Mk7PT"),
    MetVersion("UMG_Mk8PT", "2014-07-15", "2015-08-24", "MetDefnUMG_Mk8_L59PTpp.txt",
               "Use_UMG_Mk8_L59PTpp.txt", "Global/UMG_Mk8PT"),
    MetVersion("UMG_Mk9PT", "2015-08-25", "2017-07-10", "MetDefnUMG_Mk9_L59PTpp.txt",
               "Use_UMG_Mk9_L59PTpp.txt", "Global/UMG_Mk9PT"),
    MetVersion("UMG_Mk10PT", "2017-07-11", "2022-05-03", "MetDefnUMG_Mk10_L59PTpp.txt",
               "Use_UMG_Mk10_L59PTpp.txt", "Global/UMG_Mk10PT"),
    MetVersion("UMG_Mk11PT", "2022-05-04", "2026-01-20", "MetDefnUMG_Mk11_L59PTpp.txt",
               "Use_UMG_Mk11_L59PTpp.txt", "Global/UMG_Mk11PT"),
    MetVersion("UMG_Mk12PT", "2026-01-21", None, "MetDefnUMG_Mk12_L59PTpp.txt",
               "Use_UMG_Mk12_L59PTpp.txt", "Global/UMG_Mk12PT"),
)  # fmt: skip


def met_version_for_period(start, end):
    """The single met version covering every day from ``start`` to ``end`` (inclusive).

    Raises ValueError if the period spans a changeover (the run must then be shortened or split;
    not yet supported).
    """
    start = pd.Timestamp(start).normalize()
    end = pd.Timestamp(end).normalize()
    for version in MET_VERSIONS:
        first = pd.Timestamp(version.first_day)
        last = pd.Timestamp(version.last_day) if version.last_day else pd.Timestamp.max
        if first <= start and end <= last:
            return version
    raise ValueError(f"No single met version covers {start.date()} to {end.date()}")


def read_met_file_patterns(definition_path):
    """File-name pattern of every met "partition" (PT1..PT14) in a NAME met definition file.

    Returns:
        DataFrame with columns ``met_name`` (e.g. "UMG_Mk6_L59PT1pp"), ``dt_hours``, ``prefix``
        and ``suffixes`` (list of file suffixes; from Mk7 there are two files per time, "_I_" and
        "_M_").
        A met file is called <prefix><YYYYMMDDHHMM>.<suffix>.
    """
    lines = Path(definition_path).read_text().splitlines()
    header_index = next(
        i for i, line in enumerate(lines) if re.match(r"\s*Name,\s*Binary Format,", line)
    )
    headers = [h.strip() for h in lines[header_index].split(",")]
    rows = []
    for line in lines[header_index + 1 :]:
        if not line.strip():
            break
        cells = dict(zip(headers, (c.strip() for c in line.split(",")), strict=False))
        hours, minutes = cells["Dt"].split(":")
        rows.append(
            {
                "met_name": cells["Name"],
                "dt_hours": int(hours) + int(minutes) / 60,
                "prefix": cells["Prefix"],
                "suffixes": [s.strip() for s in cells["Suffix"].split(";") if s.strip()],
            }
        )
    return pd.DataFrame(rows)


def required_met_files(version, start, end, name_dir=NAME_DIR):
    """Every met file NAME needs for times from ``start`` to ``end``.

    One extra time step is included at each end, because NAME interpolates in time between met
    files on either side.

    Returns:
        Sorted list of file names (e.g. "MO201212310000.UMG_Mk6_L59PT7.pp").
    """
    patterns = read_met_file_patterns(
        Path(name_dir) / "Resources" / "Defns" / version.definition_file
    )
    dt_hours = patterns["dt_hours"].unique()
    if len(dt_hours) != 1:
        raise ValueError("Met partitions have different time steps")
    step = pd.Timedelta(hours=float(dt_hours[0]))
    first = pd.Timestamp(start).floor(step) - step
    last = pd.Timestamp(end).ceil(step) + step
    names = []
    for time in pd.date_range(first, last, freq=step):
        for _, partition in patterns.iterrows():
            for suffix in partition["suffixes"]:
                names.append(f"{partition['prefix']}{time:%Y%m%d%H%M}.{suffix}")
    return sorted(names)


def missing_met_files(version, file_names, archive=MET_ARCHIVE):
    """File names (from ``required_met_files``) found in the archive neither plain nor gzipped."""
    folder = Path(archive) / version.archive_folder
    return [
        name
        for name in file_names
        if not (folder / name).exists() and not (folder / f"{name}.gz").exists()
    ]


def _gzip_uncompressed_size(path):
    """Uncompressed size recorded in a gzip file's last 4 bytes (modulo 2**32)."""
    with open(path, "rb") as f:
        f.seek(-4, 2)
        return int.from_bytes(f.read(4), "little")


def restore_met_files(version, file_names, met_folder, archive=MET_ARCHIVE):
    """Make every met file available in ``met_folder`` before any NAME run starts.

    This replaces NAME's on-demand restore for batch runs: several NAME jobs restoring into the
    same folder at once read half-written files and lost particles (docs/jasmin.md). Here:

    * a gzipped archive file is decompressed to "<name>.partial" and then renamed, so a file
      under its final name is always complete;
    * an uncompressed archive file is linked, not copied (as the NAME team's restore script);
    * a file already present is kept if its size matches the archive, otherwise an error is
      raised (it is never silently overwritten).

    Returns:
        dict with counts: ``restored``, ``linked``, ``already_present``.
    """
    import gzip
    import os
    import shutil

    folder = Path(met_folder)
    folder.mkdir(parents=True, exist_ok=True)
    source_folder = Path(archive) / version.archive_folder
    counts = {"restored": 0, "linked": 0, "already_present": 0}
    for name in file_names:
        target = folder / name
        plain = source_folder / name
        zipped = source_folder / f"{name}.gz"
        if plain.exists():
            expected = plain.stat().st_size
        elif zipped.exists():
            expected = _gzip_uncompressed_size(zipped)
        else:
            raise FileNotFoundError(f"{name} is not in {source_folder}")

        if target.exists():
            if target.stat().st_size % 2**32 != expected % 2**32:
                raise ValueError(f"{target} exists but its size differs from the archive")
            counts["already_present"] += 1
        elif plain.exists():
            target.symlink_to(plain)
            counts["linked"] += 1
        else:
            partial = folder / f"{name}.partial"
            with gzip.open(zipped, "rb") as source, open(partial, "wb") as destination:
                shutil.copyfileobj(source, destination, length=16 * 2**20)
            if partial.stat().st_size % 2**32 != expected % 2**32:
                raise ValueError(f"Decompressed {name} has the wrong size")
            os.replace(partial, target)
            counts["restored"] += 1
    return counts
