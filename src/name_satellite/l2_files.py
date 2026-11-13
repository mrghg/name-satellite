"""Where the satellite Level-2 files are, and downloading them (input to workflow step 1).

Where this sits (see docs/workflow.md)
--------------------------------------
    CEDA  ->  **download one day's L2 file**  ->  1. ``observations.read_l2``

The University of Leicester OCPR proxy XCH4 products are open data on CEDA. The JASMIN mount
``/neodc/gosat`` is not readable for us (permission denied, group ``byacl``), so files are fetched
over HTTPS from ``dap.ceda.ac.uk`` (no login needed) into a cache folder, once.

    path = download_l2("GOSAT-1 v10.0", "2013-01-01", cache_folder)

``prepare`` uses the cache ``prepare.L2_CACHE`` = ``$NAME_SATELLITE_SCRATCH/l2_cache`` (default
``/work/scratch-pw5/$USER/name-satellite/l2_cache``). Each file is one day's netCDF as
published, e.g.
``UoL-GHG-L2-CH4-GOSAT-OCPR-20130101-fv10.0.nc``: dimensions ``n`` (soundings) × ``m`` (20
retrieval levels); ``ncdump -h`` shows its variables and the product's bias-correction note.
Nothing in the file is changed here.

GOSAT-1 and GOSAT-2 are separate products and are never mixed (decision D17).
"""

import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

CEDA = "https://dap.ceda.ac.uk/neodc/gosat/data/ch4"


@dataclass(frozen=True)
class L2Product:
    """One L2 product: its URL and file name patterns (``{date}`` = YYYYMMDD) and coverage."""

    name: str
    folder: str  # under CEDA, with {year}
    file_name: str  # with {date}
    first_day: str
    last_day: str


PRODUCTS = {
    p.name: p
    for p in (
        L2Product(
            name="GOSAT-1 v10.0",
            folder="UoL_gosat1_v10/CH4_GOSAT_OCPR/v10.0/{year}",
            file_name="UoL-GHG-L2-CH4-GOSAT-OCPR-{date}-fv10.0.nc",
            first_day="2009-04-23",
            last_day="2025-12-18",
        ),
        L2Product(
            name="GOSAT-2 v2.0",
            folder="UoL_gosat2_v2/CH4_GOSAT2_OCPR/v2.0/{year}",
            file_name="UoL-GHG-L2-CH4-GOSAT2-OCPR-{date}-fv2.0.nc",
            first_day="2019-02-05",
            last_day="2025-12-31",
        ),
    )
}


class L2FileMissing(FileNotFoundError):
    """CEDA has no file for this product and day (e.g. no soundings that day)."""


def l2_url(product_name, day):
    """URL of one day's file (``product_name``: a key of PRODUCTS; ``day``: a date).

    Raises:
        KeyError for an unknown product; ValueError outside the product's coverage.
    """
    product = PRODUCTS[product_name]
    day = pd.Timestamp(day)
    if not pd.Timestamp(product.first_day) <= day <= pd.Timestamp(product.last_day):
        raise ValueError(
            f"{day:%Y-%m-%d} is outside {product_name} ({product.first_day} to {product.last_day})"
        )
    folder = product.folder.format(year=day.year)
    return f"{CEDA}/{folder}/{product.file_name.format(date=f'{day:%Y%m%d}')}"


def download_l2(product_name, day, cache_folder):
    """Path of one day's L2 file in ``cache_folder``, downloading it first if not already there.

    The file is written to ``<name>.partial`` and renamed when complete, so an interrupted
    download never leaves a file that looks finished.

    Args:
        product_name: a key of PRODUCTS, e.g. "GOSAT-1 v10.0".
        day: the day (anything ``pandas.Timestamp`` reads).
        cache_folder: folder for downloaded files (created if needed).

    Returns:
        Path of the file, named as on CEDA.

    Raises:
        L2FileMissing if CEDA has no file for that day (HTTP 404); other HTTP errors as raised.
    """
    url = l2_url(product_name, day)
    path = Path(cache_folder) / url.rsplit("/", 1)[1]
    if path.exists():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    try:
        with urllib.request.urlopen(url, timeout=120) as response:
            partial.write_bytes(response.read())
    except urllib.error.HTTPError as error:
        if error.code == 404:
            raise L2FileMissing(
                f"No {product_name} file for {pd.Timestamp(day):%Y-%m-%d}"
            ) from error
        raise
    partial.rename(path)
    return path
