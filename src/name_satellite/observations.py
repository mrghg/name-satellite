"""Step 1 of the workflow: satellite Level-2 file → table of soundings → observation file.

Where this sits (see docs/workflow.md)
--------------------------------------
    **1. read L2, filter, (optionally) bin**  ->  2. release plan (releases.py)  ->  3. NAME input

Supported products (University of Leicester OCPR proxy XCH4; all share one file layout, dimensions
``n`` soundings × ``m`` = 20 retrieval levels, surface first):

* GOSAT-1 v10.0 and GOSAT-2 v2.0 (NCEO, on CEDA): the products for new work;
* ESA CCI GOSAT OCPR v7.x: used only to validate against legacy runs.

What comes out
--------------
``read_l2`` returns an ``xarray.Dataset`` with dimensions (``sounding``, ``level``) and only the
variables the pipeline uses, under explicit names and units:

    sounding_id (str)        the key of each entry, used to match records everywhere (never row
                             order or time; decision D18): the exposure ID of the sounding, or for
                             a legacy bin the exposure ID of its earliest sounding
    exposure_id (str)        the exposure ID(s): one, or for a bin all of them joined by commas
    time                     measurement time (see the time-resolution note below)
    latitude, longitude      degrees
    xch4, xch4_uncertainty   ppb
    xch4_quality_flag        0 = good
    retr_flag                0 = land, 1 = sun glint
    pressure_levels          hPa, (sounding, level), level 0 = surface
    pressure_weight          (sounding, level), sums to 1 over the 20 levels
    xch4_averaging_kernel    (sounding, level), normalised column averaging kernel
    ch4_profile_apriori      ppb, (sounding, level)

Filters (``quality_mask``, ``bad_pressure_mask``, ``mode_mask``, ``latlon_mask``) return boolean
masks; apply them with ``select``. ``legacy_bin`` averages soundings within lat/lon bins exactly as
the legacy ``acrg.satellite.gosat.binned_mean``.

Two writers: ``write_obs_file`` (our observation file, keyed by ``sounding_id``) and
``write_legacy_obs_file`` (the legacy layout keyed by time, for comparisons with legacy files; it
refuses repeated times).

Time resolution (docs/deviations.md X22)
----------------------------------------
GOSAT-1 v10.0 and GOSAT-2 v2.0 store the measurement time rounded DOWN to the whole minute (checked
against v7.0 for 2013-01-01: every one of the 2,024 common soundings equals the v7.0 time floored).
So many soundings share a time, and times are up to 59 s early. ``read_l2`` records this in the
``time_resolution`` attribute. The legacy preprocessing cannot handle these files at all (its
binning needs unique times). Our files are keyed by sounding ID instead (D18).

Legacy source: ``acrg/satellite/gosat.py`` (filters ~103-270, binning ~425-886,
``gosat_process_file`` ~1191) and ``common.py`` (``add_coords``, ``output``).
"""

from pathlib import Path

import numpy as np
import xarray as xr

from name_satellite.releases import footprint_box_degrees

N_LEVELS = 20
BAD_PRESSURE = -9999.99  # legacy fill value for an undefined pressure level

# L2 variable -> our name (most are unchanged; the L2 files call the weights "pressure_weight").
L2_VARIABLES = {
    "exposure_id": "exposure_id",
    "time": "time",
    "latitude": "latitude",
    "longitude": "longitude",
    "xch4": "xch4",
    "xch4_uncertainty": "xch4_uncertainty",
    "xch4_quality_flag": "xch4_quality_flag",
    "retr_flag": "retr_flag",
    "pressure_levels": "pressure_levels",
    "pressure_weight": "pressure_weight",
    "xch4_averaging_kernel": "xch4_averaging_kernel",
    "ch4_profile_apriori": "ch4_profile_apriori",
}
PROFILE_VARIABLES = (
    "pressure_levels",
    "pressure_weight",
    "xch4_averaging_kernel",
    "ch4_profile_apriori",
)


# ---------------------------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------------------------


def identify_product(attrs, filename):
    """Return (product, version) of a UoL OCPR file from its global attributes and file name.

    Raises ValueError for anything else, rather than guessing.
    """
    title = attrs.get("title", "")
    if "ESA CCI GOSAT OCPR CH4" in title:
        return "ESA CCI GOSAT OCPR", attrs.get("product_version", "")
    if "GOSAT2" in filename and "OCPR" in filename:
        return "UoL GOSAT-2 OCPR", filename.rsplit("-fv", 1)[-1].removesuffix(".nc")
    if "GOSAT" in filename and "OCPR" in filename:
        return "UoL GOSAT OCPR", filename.rsplit("-fv", 1)[-1].removesuffix(".nc")
    raise ValueError(f"Not a recognised UoL OCPR proxy XCH4 file: {filename} (title {title!r})")


def read_l2(path):
    """Read one day's UoL OCPR proxy XCH4 file into a sounding table.

    Args:
        path: L2 netCDF file (GOSAT-1 v10, GOSAT-2 v2, or ESA CCI v7.x).

    Returns:
        xarray.Dataset with dims (sounding, level); variables as in the module docstring, in the
        file's own units (ppb, hPa). Soundings keep the file's order. Global attributes record the
        product, version, source file, the bias correction note (if any) and ``time_resolution``.

    Raises:
        ValueError if the layout is not the expected (n soundings × 20 levels), exposure IDs are not
        unique, or a level's pressure is not below the one beneath it.
    """
    path = Path(path)
    raw = xr.open_dataset(path).load()
    product, version = identify_product(raw.attrs, path.name)
    if raw.sizes.get("m") != N_LEVELS or "n" not in raw.sizes:
        raise ValueError(f"{path.name}: expected dimensions n × m=20, got {dict(raw.sizes)}")

    data = {}
    for l2_name, name in L2_VARIABLES.items():
        values = raw[l2_name].values
        if l2_name in PROFILE_VARIABLES:
            data[name] = (("sounding", "level"), values)
        else:
            data[name] = ("sounding", values)
    table = xr.Dataset(data)
    table["exposure_id"] = ("sounding", [x.decode().strip() for x in raw["exposure_id"].values])
    table["sounding_id"] = table["exposure_id"].copy()  # D18: the sounding's ID is its key

    if len(set(table["exposure_id"].values)) != table.sizes["sounding"]:
        raise ValueError(f"{path.name}: exposure IDs are not unique")
    # Level 0 must be the surface: if most soundings were ordered top first, the file is not laid
    # out as we assume. (A few soundings over high ground have level 0 above level 1: see
    # level_order_mask.)
    pressure = table["pressure_levels"].values
    if np.nanmedian(pressure[:, 0]) < np.nanmedian(pressure[:, -1]):
        raise ValueError(f"{path.name}: pressure levels are not ordered surface first")

    times = table["time"].to_index()
    whole_minutes = bool((times.second == 0).all() and (times.microsecond == 0).all())
    table.attrs = {
        "product": product,
        "product_version": version,
        "source_file": path.name,
        "bias_correction": raw.attrs.get("bias_correction", "none stated"),
        "time_resolution": (
            "whole minutes, rounded down (docs/deviations.md X22)"
            if whole_minutes
            else "seconds or finer"
        ),
    }
    return table


# ---------------------------------------------------------------------------------------------
# Filters: each returns a boolean mask over soundings (True = keep)
# ---------------------------------------------------------------------------------------------


def quality_mask(table):
    """Keep soundings flagged good (``xch4_quality_flag == 0``). Legacy ``gosat_quality_filter``."""
    return table["xch4_quality_flag"].values == 0


def _pressure_defined(pressure):
    """True for soundings whose 20 pressure levels are all defined (not NaN, not -9999.99)."""
    return ~(np.isnan(pressure) | (pressure == np.float32(BAD_PRESSURE))).any(axis=1)


def bad_pressure_mask(table):
    """Drop soundings with any undefined pressure level.

    Legacy ``gosat_pressure_filter`` looks for -9999.99 only. In files that declare -9999.99 as the
    netCDF fill value (GOSAT-2 v2.0), xarray reads those levels as NaN, so the legacy filter removes
    nothing (docs/deviations.md X24). We treat NaN as undefined too. In files where this filter
    also matched in legacy, the result is the same.
    """
    return _pressure_defined(table["pressure_levels"].values)


def level_order_mask(table):
    """Keep soundings whose pressure decreases strictly from level to level (surface first).

    Over high ground (e.g. Tibet, Antarctica) the retrieved surface pressure (level 0) can be LOWER
    than level 1, which would give the bottom release layer a negative thickness. Legacy removed
    such soundings by comparing level 1 with a NAME surface pressure ("level_order" in
    ``name_pressure_filter``); this uses the retrieval's own levels. Soundings with undefined levels
    are dropped too.
    """
    pressure = table["pressure_levels"].values
    with np.errstate(invalid="ignore"):
        decreasing = (np.diff(pressure, axis=1) < 0).all(axis=1)
    return decreasing & _pressure_defined(pressure)


def mode_mask(table, mode):
    """Keep land (``retr_flag == 0``) or sun-glint (``retr_flag == 1``) soundings."""
    flags = {"land": 0, "glint": 1}
    if mode not in flags:
        raise ValueError(f"mode must be 'land' or 'glint', not {mode!r}")
    return table["retr_flag"].values == flags[mode]


def latlon_mask(table, lat_bounds, lon_bounds):
    """Keep soundings with lat_bounds[0] <= latitude < lat_bounds[1] (same for longitude).

    Lower bound included, upper excluded, as legacy ``latlon_filter``.
    """
    lat = table["latitude"].values
    lon = table["longitude"].values
    return (
        (lat >= lat_bounds[0])
        & (lat < lat_bounds[1])
        & (lon >= lon_bounds[0])
        & (lon < lon_bounds[1])
    )


def select(table, mask):
    """The soundings where ``mask`` is True (attributes kept)."""
    return table.isel(sounding=np.flatnonzero(mask))


# ---------------------------------------------------------------------------------------------
# Legacy binning
# ---------------------------------------------------------------------------------------------


def legacy_bin_indices(table, lat_bounds, lon_bounds, bin_size):
    """Bin index (lat, lon) of each sounding, or -1 for soundings outside the bins.

    Legacy ``create_zbin``: bin edges ``np.arange(lower, upper + size, size)``, positions assigned
    with ``np.digitize``; soundings below the first edge or at/after the last are dropped
    (``zbin_filter``).

    Args:
        bin_size: (lat, lon) bin size in degrees, e.g. (0.234, 0.352).

    Returns:
        (lat_index, lon_index) integer arrays; -1 where the sounding is outside.
    """
    indices = []
    for values, bounds, size in (
        (table["latitude"].values, lat_bounds, bin_size[0]),
        (table["longitude"].values, lon_bounds, bin_size[1]),
    ):
        edges = np.arange(bounds[0], bounds[1] + size, size)
        index = np.digitize(values, edges)  # 0 = below the first edge, len(edges) = above the last
        index = np.where((index == 0) | (index == len(edges)), -1, index - 1)
        indices.append(index)
    return indices[0], indices[1]


def legacy_bin(table, lat_bounds, lon_bounds, bin_size, footprint_km=10.5):
    """Average soundings that fall in the same lat/lon bin, as legacy ``binned_mean``.

    Within each bin (legacy ``mean_ds``; LEGACY QUIRKS, see docs/deviations.md X23):
      * every variable is the plain mean over the bin's soundings (no weighting);
      * time is the midpoint of the first and last sounding;
      * ``xch4_uncertainty`` = sqrt(mean(uncertainty)² + std(xch4)²) (population std): it does not
        decrease with the number of soundings averaged;
      * the release box is the footprint size plus the spread of the positions:
        dlat = footprint + (max - min latitude), dlon = footprint at the mean latitude
        + (max - min longitude);
      * exposure IDs are joined with commas; the bin's ``sounding_id`` is the exposure ID of its
        earliest sounding (D18, not legacy).
    The legacy code also adds random microsecond offsets to duplicate times (X12); this does not.

    Args:
        table: sounding table (after filtering).
        lat_bounds, lon_bounds: (lower, upper) degrees.
        bin_size: (lat, lon) bin size in degrees.
        footprint_km: sounding footprint diameter (km); 10.5 for GOSAT-1.

    Returns:
        Dataset like ``table``, one entry per occupied bin, sorted by time, plus ``dlat`` and
        ``dlon`` (degrees, full widths) and ``n_soundings``.
    """
    lat_index, lon_index = legacy_bin_indices(table, lat_bounds, lon_bounds, bin_size)
    inside = (lat_index >= 0) & (lon_index >= 0)
    keys = sorted(set(zip(lat_index[inside], lon_index[inside], strict=True)))

    rows = []
    for key in keys:
        members = np.flatnonzero(inside & (lat_index == key[0]) & (lon_index == key[1]))
        rows.append(_bin_one(table.isel(sounding=members), footprint_km))
    binned = xr.concat(rows, dim="sounding")
    binned = binned.isel(sounding=np.argsort(binned["time"].values, kind="stable"))
    binned.attrs = dict(table.attrs)
    binned.attrs["binning"] = (
        f"legacy_bin: latitude {lat_bounds[0]} to {lat_bounds[1]}, longitude {lon_bounds[0]} to "
        f"{lon_bounds[1]}, bins of {bin_size[0]} x {bin_size[1]} degrees (lat x lon)"
    )
    return binned


def _bin_one(members, footprint_km):
    """Legacy mean of one bin's soundings (see ``legacy_bin``)."""
    out = {}
    # D18: a bin is identified by the exposure ID of its earliest sounding (ties: smallest ID).
    earliest = np.lexsort((members["exposure_id"].values, members["time"].values))[0]
    for name in members.data_vars:
        values = members[name]
        if name == "sounding_id":
            out[name] = str(members["exposure_id"].values[earliest])
        elif name == "exposure_id":
            out[name] = ",".join(values.values)
        elif name == "time":
            times = np.sort(values.values)
            out[name] = times[-1] - (times[-1] - times[0]) / 2
        elif name == "xch4_uncertainty":
            out[name] = np.sqrt(np.mean(values.values) ** 2 + np.std(members["xch4"].values) ** 2)
        else:
            out[name] = values.mean("sounding").values
    lat = members["latitude"].values
    lon = members["longitude"].values
    dlon_footprint, dlat_footprint = footprint_box_degrees(np.mean(lat), footprint_km)
    out["dlat"] = dlat_footprint + (lat.max() - lat.min())
    out["dlon"] = dlon_footprint + (lon.max() - lon.min())
    out["n_soundings"] = members.sizes["sounding"]

    data = {}
    for name, value in out.items():
        if name in PROFILE_VARIABLES:
            data[name] = (("sounding", "level"), np.asarray(value)[None, :])
        else:
            data[name] = ("sounding", [value])
    return xr.Dataset(data)


# ---------------------------------------------------------------------------------------------
# Legacy observation file (read by the inversion code and by releases.build_release_plan)
# ---------------------------------------------------------------------------------------------

# our name -> legacy file name (acrg.satellite.common.output)
LEGACY_OBS_NAMES = {
    "xch4": "xch4",
    "xch4_uncertainty": "xch4_uncertainty",
    "latitude": "lat",
    "longitude": "lon",
    "pressure_levels": "pressure_levels",
    "pressure_weight": "pressure_weights",
    "xch4_averaging_kernel": "xch4_averaging_kernel",
    "ch4_profile_apriori": "ch4_profile_apriori",
}


def to_legacy_obs(table):
    """The sounding table in the legacy observation-file layout (one day).

    Dims ``time`` and ``lev`` (0-19), plus ``id`` for the exposure IDs: ``exposure_id`` is
    (time, id), one column per sounding combined into a bin, padded with "nan" as legacy
    (``split_output``). Variable names as legacy (``lat``, ``lon``, ``pressure_weights``).
    ``dlat``/``dlon`` of binned soundings are NOT part of the legacy file; pass them to
    ``releases.build_release_plan`` as ``box_degrees`` instead.

    Raises:
        ValueError if two soundings have the same time: the legacy file and the inversion code
        index by time (docs/deviations.md X22, X12).
    """
    times = table["time"].values
    if len(np.unique(times)) != len(times):
        raise ValueError(
            "Duplicate sounding times: the legacy observation file is indexed by time "
            "(see docs/deviations.md X22)"
        )
    data = {}
    for name, legacy_name in LEGACY_OBS_NAMES.items():
        if name in PROFILE_VARIABLES:
            data[legacy_name] = (("time", "lev"), table[name].values)
        else:
            data[legacy_name] = ("time", table[name].values)
    ids = [value.split(",") for value in table["exposure_id"].values]
    width = max(len(i) for i in ids)
    padded = np.array([i + ["nan"] * (width - len(i)) for i in ids], dtype="S")
    data["exposure_id"] = (("time", "id"), padded)
    coords = {"time": times, "lev": np.arange(N_LEVELS), "id": np.arange(1, width + 1)}
    obs = xr.Dataset(data, coords=coords, attrs=dict(table.attrs))
    for name, legacy_name in LEGACY_OBS_NAMES.items():
        obs[legacy_name].attrs["units"] = _UNITS[name]
    return obs


_UNITS = {
    "xch4": "1e-9",
    "xch4_uncertainty": "1e-9",
    "latitude": "degrees_north",
    "longitude": "degrees_east",
    "pressure_levels": "hPa",
    "pressure_weight": "1",
    "xch4_averaging_kernel": "1",
    "ch4_profile_apriori": "1e-9",
}


def write_legacy_obs_file(table, path):
    """Write one day's soundings as a legacy observation file (``to_legacy_obs`` layout).

    The legacy file name is ``gosat-fts_gosat_YYYYMMDD_ch4-column.nc``.
    """
    to_legacy_obs(table).to_netcdf(path)


# ---------------------------------------------------------------------------------------------
# Our observation file: keyed by sounding ID (decision D18)
# ---------------------------------------------------------------------------------------------


def to_obs(table):
    """The sounding table as our observation file (one day, one satellite product: D17).

    Like the legacy layout (``to_legacy_obs``: same variable names and units, ``lev`` 0-19, and
    ``exposure_id`` as (sounding, id) padded with "nan"), but the soundings are indexed by
    ``sounding_id`` instead of time, and ``time`` is an ordinary variable. GOSAT-1 v10 and GOSAT-2
    times repeat (whole minutes, X22), so time cannot be the key. Binned soundings also carry
    ``dlat``/``dlon`` (release box, degrees) and ``n_soundings``.

    Raises:
        ValueError if a sounding ID is repeated.
    """
    ids = table["sounding_id"].values
    if len(set(ids)) != len(ids):
        raise ValueError("Repeated sounding IDs: they must identify the soundings uniquely")
    dims = "sounding_id"
    data = {"time": (dims, table["time"].values)}
    for name, legacy_name in LEGACY_OBS_NAMES.items():
        if name in PROFILE_VARIABLES:
            data[legacy_name] = ((dims, "lev"), table[name].values)
        else:
            data[legacy_name] = (dims, table[name].values)
    for name in ("dlat", "dlon", "n_soundings"):
        if name in table:
            data[name] = (dims, table[name].values)
    exposures = [value.split(",") for value in table["exposure_id"].values]
    width = max(len(i) for i in exposures)
    padded = np.array([i + ["nan"] * (width - len(i)) for i in exposures], dtype="S")
    data["exposure_id"] = ((dims, "id"), padded)
    coords = {
        dims: np.array(ids, dtype=str),
        "lev": np.arange(N_LEVELS),
        "id": np.arange(1, width + 1),
    }
    obs = xr.Dataset(data, coords=coords, attrs=dict(table.attrs))
    for name, legacy_name in LEGACY_OBS_NAMES.items():
        obs[legacy_name].attrs["units"] = _UNITS[name]
    for name in ("dlat", "dlon"):
        if name in obs:
            obs[name].attrs["units"] = "degrees"
            obs[name].attrs["long_name"] = f"release box full width ({name[1:]}itude)"
    obs["sounding_id"].attrs["long_name"] = (
        "sounding ID: exposure ID of the sounding (or of the earliest sounding of a bin)"
    )
    return obs


def write_obs_file(table, path):
    """Write one day's soundings as our observation file (``to_obs`` layout)."""
    to_obs(table).to_netcdf(path)
