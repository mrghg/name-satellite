"""Registry of legacy reference data used to validate the pipeline (LEGACY VALIDATION ONLY).

Not part of the workflow for new runs; used by the tests and validation scripts.

The files themselves are NOT in git (they live under ``data/``, or wherever the
``NAME_SATELLITE_REFDATA`` environment variable points). This module records
where each file is expected and its SHA-256 checksum, so tests can check that
they are validating against exactly the files described in
docs/reference-data.md.
"""

from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass
from pathlib import Path

ENV_VAR = "NAME_SATELLITE_REFDATA"

# Run folders, relative to the reference-data root.
# Full January 2013 of the 2018 concentration run (supersedes the earlier 1-day copy).
CONC_2018 = "gosat_brazil_201301/SOUTHAMERICA_GOSAT-BRAZIL_column_conc"
CONC_2018_MANIFEST = "dev/manifests/gosat_brazil_201301_conc.sha256"
MR_2019 = "gosat_brazil_20130101_mr_2019/SOUTHAMERICA_GOSAT-BRAZIL_column"
# Full January 2013 of the 2019 mixing-ratio run (the production run). Its 1 Jan files are the same
# as MR_2019's; MR_2019 is kept as it also holds time_step.txt and the NAME_processing benchmark.
MR_2019_MONTH = "gosat_brazil_201301_mr_2019/SOUTHAMERICA_GOSAT-BRAZIL_column"
MR_2019_MONTH_MANIFEST = "dev/manifests/gosat_brazil_201301_mr_2019.sha256"
RUN = "SOUTHAMERICA_GOSAT-BRAZIL_column_20130101"
# ESA CCI v7.0 L2 files and the legacy acrg.satellite preprocessing of them (observations tests).
CCI_V70 = "gosat_cci_v7.0_201301_legacy_obs"


@dataclass(frozen=True)
class RefFile:
    key: str
    relpath: str
    sha256: str
    description: str


REFERENCE_FILES: tuple[RefFile, ...] = (
    # --- April 2019 mixing-ratio run (the one production fp was made from) ---
    RefFile(
        "mr2019_input",
        f"{MR_2019}/Input_files/BackRun_{RUN}.txt",
        "ff1eb3eecbda3d9b9b7b4e577a862d6472a0b0b4f9227b2cad452ce7508d4f94",
        "NAME 7.2 main input file, 2019 mixing-ratio run",
    ),
    RefFile(
        "mr2019_met_input",
        f"{MR_2019}/Input_files/BackRun_{RUN}_met.txt",
        "e5b24aabaf1f0e5e43584845aa8986458d5c5239baeb92e9e1446d6381e8de00",
        "NAME 7.2 met-extraction input file, 2019 run",
    ),
    RefFile(
        "mr2019_fields",
        f"{MR_2019}/MixR_files/Fields_{RUN}.txt.gz",
        "74864f6957c28c9c996f80061deb670137c36e9d0836806f185350065f5b3979",
        "Time-integrated mixing-ratio fields, 425 field columns",
    ),
    RefFile(
        "mr2019_particles",
        f"{MR_2019}/Particle_files/particle_location_{RUN}.txt.gz",
        "01c4e8d0031fb05014d98f2659bc20f038fc2fc25589d4c2943b904d9e3a671c",
        "Particle fate file (ACRG format), 562,925 particles",
    ),
    RefFile(
        "obs",
        f"{MR_2019}/Observations/gosat-fts_gosat_20130101_ch4-column.nc",
        "ab1a077e9ee3547d6b67ca725532dd56d44abebde3e2a1dff1e8c5c2f2f846dd",
        "Binned GOSAT OCPR v7.2 obs file (25 points); identical in both runs",
    ),
    RefFile(
        "mr2019_benchmark",
        f"{MR_2019}/Processed_Fields_files_Benchmark/GOSAT-BRAZIL-column_SOUTHAMERICA_201301.nc",
        "d73458e866448ca6ff0917945ecd99041f370a6e9081c36b9412fa88630f6e55",
        "Processed fp for 2013-01-01 only (from NAME_processing tests)",
    ),
    RefFile(
        "mr2019_month_fields_20130131",
        f"{MR_2019_MONTH}/MixR_files/Fields_SOUTHAMERICA_GOSAT-BRAZIL_column_20130131.txt.gz",
        "ef6ed0af65818595a6520f6ed689c6ed8a011f1b7407e2413c6bbe3409b3f451",
        "Mixing-ratio fields of the last day of the month copy (full manifest in dev/manifests)",
    ),
    # --- ESA CCI v7.0 L2 files and their legacy preprocessing (observations.py tests) ---
    RefFile(
        "cci_v70_l2_20130101",
        f"{CCI_V70}/l2/ESACCI-GHG-L2-CH4-GOSAT-OCPR-20130101-fv7.0.nc",
        "a0e686c34ec48c0049868104f3c847a924149a5d1a64335682fe8d2360796ff6",
        "ESA CCI GOSAT OCPR v7.0 L2 file, 2013-01-01 (from /neodc/esacci)",
    ),
    RefFile(
        "cci_v70_l2_20130107",
        f"{CCI_V70}/l2/ESACCI-GHG-L2-CH4-GOSAT-OCPR-20130107-fv7.0.nc",
        "7e5722f8d862ac4c4b29c52520cda7ea1821006e11606ec62a471e79ac4a8251",
        "ESA CCI GOSAT OCPR v7.0 L2 file, 2013-01-07",
    ),
    RefFile(
        "cci_v70_legacy_binned_20130101",
        f"{CCI_V70}/legacy_20130101/binned.nc",
        "79f06adc5b8061fb3a46de1667a791794ce15ce26985dfcf8414401c9bf5424d",
        "Legacy gosat_process_file output, binned (legacy_validation/run_legacy_gosat.py)",
    ),
    RefFile(
        "cci_v70_legacy_unbinned_20130101",
        f"{CCI_V70}/legacy_20130101/unbinned.nc",
        "0f70999b0ff09a06350c70c67155c412d14214ba20ac4eee199a97370e86ac9d",
        "Legacy gosat_process_file output, lat/lon filter only",
    ),
    RefFile(
        "cci_v70_legacy_binned_20130107",
        f"{CCI_V70}/legacy_20130107/binned.nc",
        "ed3b0f685ec9f8fa3ef2322e69c70442dd7c9fa79a3d70e0869d4a551fff1068",
        "Legacy gosat_process_file output, binned",
    ),
    RefFile(
        "cci_v70_legacy_unbinned_20130107",
        f"{CCI_V70}/legacy_20130107/unbinned.nc",
        "9b491e3bec4a91e61e19df8e5c1a9e5bb451db4c11b539f932e9ea8cf27b0c4b",
        "Legacy gosat_process_file output, lat/lon filter only",
    ),
    # --- Production processed file (whole month) ---
    RefFile(
        "production_201301",
        "GOSAT-BRAZIL-column_SOUTHAMERICA_201301.nc",
        "ac92148e3d8209864b2ccadb2635df44d0171015a330381598201ba1b2ca45ae",
        "Production fp file, Jan 2013 (678 obs); 2013-01-01 slice == benchmark",
    ),
    # --- June 2018 concentration run (superseded; for noise/unit comparisons) ---
    RefFile(
        "conc2018_input",
        f"{CONC_2018}/Input_files/BackRun_{RUN}.txt",
        "a3a245b2530706a50e8b3b55836fc9d7bb39a3e83b2b33e14987e1a1079f7234",
        "NAME 7.2 main input file, 2018 air-concentration run",
    ),
    RefFile(
        "conc2018_fields",
        f"{CONC_2018}/Fields_files/Fields_{RUN}.txt.gz",
        "05be0c9adba33317d154213bd2c18da31a938549d0cd6698d8cebeaae11b8233",
        "Time-integrated air-concentration fields (g s/m3)",
    ),
    RefFile(
        "conc2018_particles",
        f"{CONC_2018}/Particle_files/Original/particle_location_{RUN}.txt.gz",
        "1e4c0df81c20ef87cebf18b0451d9313b13d0f65902ed1d6d8f1ff3c5f568e98",
        "Particle fate file as written by NAME (Reprocessed = whitespace only)",
    ),
    RefFile(
        "conc2018_processed_201301",
        f"{CONC_2018}/Processed_Fields_files/GOSAT-BRAZIL-column_SOUTHAMERICA_201301.nc",
        "a4d40bc577e6813f1515be70dcd4ef0537902a8fb6b7fb093c772170a25c38e1",
        "fp processed from the conc run (created 2018-09-14); superseded, fp ~1.42x production",
    ),
)


def refdata_root() -> Path:
    """Root directory of the reference data (env var, else ``<repo>/data``)."""
    env = os.environ.get(ENV_VAR)
    if env:
        return Path(env).expanduser()
    return Path(__file__).resolve().parents[3] / "data"


def ref_path(key: str) -> Path:
    """Absolute path of a registered reference file (it may not exist)."""
    for ref in REFERENCE_FILES:
        if ref.key == key:
            return refdata_root() / ref.relpath
    raise KeyError(f"Unknown reference file key: {key!r}")


def sha256sum(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()
