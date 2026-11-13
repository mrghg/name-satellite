"""Reproduce the LEGACY satellite processing for one day of NAME output (LEGACY VALIDATION ONLY).

Where this sits (see docs/workflow.md)
--------------------------------------
This chains workflow steps 5 (read NAME output) and 6 (footprints, boundary histograms, column
combination) in exactly the way the legacy ACRG code did, so that its output can be compared with
the legacy footprint files bit for bit (docs/reference-data.md). Every legacy quirk is reproduced
on purpose and labelled with its docs/deviations.md ID. New runs are processed by
``processing.py`` instead (D13, D23, D24); this module stays as the reference.

Inputs (one legacy run folder, one day)
---------------------------------------
    <run>/MixR_files/Fields_<run>_<date>.txt.gz           time-integrated mixing ratio (ppm s)
    <run>/Particle_files/particle_location_<run>_<date>.txt.gz
    <run>/Met/Met_<run>_<date>/Met_*_<PPPLL>_C1_X1_Y1.txt.gz   release-point met
    <run>/Observations/gosat-fts_gosat_<date>_ch4-column.nc    retrieval (levels, weights, AK)
    <run>/time_step.txt                                    release length in seconds (60)

The NAME sources are ``SourceID_1 .. SourceID_N``, ordered sounding by sounding with the 17 levels
of each sounding in turn (this ordering is checked against the met files' PPPLL labels).

Output
------
An xarray.Dataset with one entry per sounding (dimension ``time``), using the legacy variable names
and units (see ``write_legacy_fp_file``).
"""

from pathlib import Path

import numpy as np
import xarray as xr

from name_satellite.boundary_histograms import legacy_boundary_histograms
from name_satellite.column import column_weights, combine_levels, nearest_release_levels
from name_satellite.footprints import footprint_from_mixing_ratio, grid_cell_area_m2
from name_satellite.legacy import legacy_grid_coordinates, read_legacy_met_folder
from name_satellite.name_output import read_acrg_particle_file, read_fields

LEVELS_PER_SOUNDING = 17  # legacy upper_level


def process_day_legacy(run_folder, run_name, date, max_level=17, histogram_variant="2019"):
    """Legacy processing of one day of satellite NAME output.

    Args:
        run_folder: the legacy run folder (contains MixR_files, Particle_files, Met, ...).
        run_name: e.g. "SOUTHAMERICA_GOSAT-BRAZIL_column".
        date: "YYYYMMDD".
        max_level: retrieval levels 0..max_level-1 are combined (legacy: 17).
        histogram_variant: "2019" (the code that made the 2019 production files) or "2025"
            (today's legacy code); see boundary_histograms.

    Returns:
        xarray.Dataset with dims (time, lat, lon, height). See ``write_legacy_fp_file`` for names.
    """
    run_folder = Path(run_folder)
    tag = f"{run_name}_{date}"

    # --- Step 5: read NAME output and the matching observations ---------------------------------
    fields = read_fields(run_folder / "MixR_files" / f"Fields_{tag}.txt.gz")
    if set(fields["units"].values) != {"ppm s"}:
        raise ValueError(f"Expected mixing ratio in ppm s, got {set(fields['units'].values)}")
    n_sources = fields.sizes["source"]
    expected_names = [f"SourceID_{n}" for n in range(1, n_sources + 1)]
    if fields["source"].values.tolist() != expected_names:
        raise ValueError("Fields columns are not SourceID_1..N in order")
    if n_sources % LEVELS_PER_SOUNDING:
        raise ValueError(f"{n_sources} sources is not a multiple of {LEVELS_PER_SOUNDING} levels")
    n_soundings = n_sources // LEVELS_PER_SOUNDING

    met = read_legacy_met_folder(run_folder / "Met" / f"Met_{tag}")
    source_from_met = (met["point_number"] - 1) * LEVELS_PER_SOUNDING + met["level_number"]
    if source_from_met.tolist() != list(range(1, n_sources + 1)):
        raise ValueError("Met files do not match the sources one-to-one in order")

    particles = read_acrg_particle_file(
        run_folder / "Particle_files" / f"particle_location_{tag}.txt.gz"
    )

    obs_files = sorted((run_folder / "Observations").glob(f"*_{date}_*.nc"))
    if len(obs_files) != 1:
        raise ValueError(f"Expected one observation file for {date}, found {obs_files}")
    obs = xr.open_dataset(obs_files[0]).load()
    if obs.sizes["time"] != n_soundings:
        raise ValueError(f"{obs.sizes['time']} soundings in obs file but {n_soundings} in NAME run")

    release_seconds = float((run_folder / "time_step.txt").read_text())
    return process_legacy_inputs(
        fields, particles, met, obs, release_seconds, max_level, histogram_variant
    )


def process_legacy_inputs(
    fields, particles, met, obs, release_seconds, max_level=17, histogram_variant="2019"
):
    """Legacy processing (workflow step 6) of one day's NAME output already read into memory.

    This is the part of ``process_day_legacy`` after the files are read. It is separate so that
    new NAME runs (e.g. the NAME 8.7.1 validation runs, whose files are named and laid out
    differently) go through exactly the same legacy processing.

    Args:
        fields: ``read_fields`` output in ppm s, sources ``SourceID_1..N`` in order.
        particles: one row per particle with ``lon``, ``lat`` (degrees), ``height_m_agl`` (m)
            and ``source_number`` (as from ``read_acrg_particle_file``).
        met: release-point met, one row per source in source order
            (``read_legacy_met_folder``).
        obs: retrieval for the same day (pressure levels, weights, averaging kernel).
        release_seconds: length of each release (s); legacy 60.
        max_level, histogram_variant: as for ``process_day_legacy``.

    Returns:
        As ``process_day_legacy``.
    """
    if set(fields["units"].values) != {"ppm s"}:
        raise ValueError(f"Expected mixing ratio in ppm s, got {set(fields['units'].values)}")
    n_sources = fields.sizes["source"]
    if fields["source"].values.tolist() != [f"SourceID_{n}" for n in range(1, n_sources + 1)]:
        raise ValueError("Fields columns are not SourceID_1..N in order")
    n_soundings = n_sources // LEVELS_PER_SOUNDING
    source_from_met = (met["point_number"] - 1) * LEVELS_PER_SOUNDING + met["level_number"]
    if source_from_met.tolist() != list(range(1, n_sources + 1)):
        raise ValueError("Met files do not match the sources one-to-one in order")
    if obs.sizes["time"] != n_soundings:
        raise ValueError(f"{obs.sizes['time']} soundings in obs file but {n_soundings} in NAME run")

    # Legacy stores the release length as hours (seconds / 3600) and multiplies back by 3600.
    release_duration_s = 3600.0 * (release_seconds / 3600.0)

    # --- Step 6a: footprints for every release ---------------------------------------------------
    lat = legacy_grid_coordinates(
        fields.attrs["Y grid origin"], fields.attrs["Y grid resolution"], fields.sizes["lat"]
    )
    lon = legacy_grid_coordinates(
        fields.attrs["X grid origin"], fields.attrs["X grid resolution"], fields.sizes["lon"]
    )
    area = grid_cell_area_m2(lat, lon)
    footprints = footprint_from_mixing_ratio(fields.values, area, release_duration_s)

    # --- Step 6b: boundary histograms for every release ------------------------------------------
    histograms = legacy_boundary_histograms(particles, lat, lon, n_sources, histogram_variant)
    edges = ("north", "east", "south", "west")

    # --- Step 6c: combine the 17 releases of each sounding into one column -----------------------
    column_fp = np.zeros((n_soundings, len(lat), len(lon)))
    column_hist = {edge: np.zeros((n_soundings, *histograms[edge].shape[1:])) for edge in edges}
    for p in range(n_soundings):
        releases = slice(p * LEVELS_PER_SOUNDING, (p + 1) * LEVELS_PER_SOUNDING)
        release_lon = met["lon"].iloc[releases.start]
        release_lat = met["lat"].iloc[releases.start]
        # Legacy only logs this; we stop, as it would mean the files are mismatched.
        if (
            abs(obs["lon"].values[p] - release_lon) > 1
            or abs(obs["lat"].values[p] - release_lat) > 1
        ):
            raise ValueError(f"Sounding {p + 1}: observation and release positions differ by >1°")

        pw = obs["pressure_weights"].values[p]
        ak = obs["xch4_averaging_kernel"].values[p]
        chosen = nearest_release_levels(
            met["pressure_pa"].to_numpy()[releases], obs["pressure_levels"].values[p], max_level
        )
        column_fp[p] = combine_levels(footprints[releases], chosen, pw, ak)
        for edge in edges:
            column_hist[edge][p] = combine_levels(histograms[edge][releases], chosen, pw, ak)

        # Legacy consistency check: combined histograms sum to Σ (pressure weight × AK).
        total = sum(column_hist[edge][p].sum() for edge in edges)
        if not np.allclose(total, column_weights(pw, ak, max_level).sum()):
            raise ValueError(f"Sounding {p + 1}: histogram total does not match Σ PW×AK")

    # --- Release-point met of each sounding: legacy keeps that of the lowest release (level 1) ---
    lowest = met.iloc[::LEVELS_PER_SOUNDING].reset_index(drop=True)

    result = xr.Dataset(
        {
            "fp": (("time", "lat", "lon"), column_fp),
            "particle_locations_n": (("time", "lon", "height"), column_hist["north"]),
            "particle_locations_e": (("time", "lat", "height"), column_hist["east"]),
            "particle_locations_s": (("time", "lon", "height"), column_hist["south"]),
            "particle_locations_w": (("time", "lat", "height"), column_hist["west"]),
            # X17: legacy writes temperature in °C but labels it K; reproduced here.
            "temperature": ("time", lowest["temperature_c"].to_numpy()),
            "pressure": ("time", lowest["pressure_pa"].to_numpy() / 100.0),  # Pa -> hPa
            "wind_speed": ("time", lowest["wind_speed_m_s"].to_numpy()),
            "wind_direction": ("time", lowest["wind_direction_deg"].to_numpy()),
            "PBLH": ("time", lowest["boundary_layer_depth_m"].to_numpy()),
            "release_lon": ("time", lowest["lon"].to_numpy()),
            "release_lat": ("time", lowest["lat"].to_numpy()),
        },
        # Legacy replaces NAME's minute-precision times with the observation times.
        coords={"time": obs["time"].values, "lat": lat, "lon": lon, "height": histograms["height"]},
        attrs={"max_level": max_level},
    )

    # Legacy drops soundings whose footprint AND histograms are all zero.
    fp_total = result["fp"].sum(("lat", "lon"))
    histogram_total = 0
    for name in HISTOGRAM_NAMES:
        other_dims = [dim for dim in result[name].dims if dim != "time"]
        histogram_total = histogram_total + result[name].sum(other_dims)
    keep = (fp_total > 0) | (histogram_total > 0)
    return result.isel(time=np.flatnonzero(keep.values))


HISTOGRAM_NAMES = (
    "particle_locations_n",
    "particle_locations_e",
    "particle_locations_s",
    "particle_locations_w",
)

# Variable layout of the legacy footprint file (acrg 75793919 ``write_netcdf``):
# name -> (netCDF dimensions, least_significant_digit, units). least_significant_digit makes the
# netCDF library round the data before compression (lossy, docs/deviations.md X7); we reproduce it
# so that our files can be compared with the legacy ones exactly.
LEGACY_FILE_VARIABLES = {
    "fp": (("lat", "lon", "time"), 5, "(mol/mol)/(mol/m2/s)"),
    "temperature": (("time",), 4, "K"),  # values are actually °C (X17)
    "pressure": (("time",), 4, "hPa"),
    "wind_speed": (("time",), 4, "m/s"),
    "wind_direction": (("time",), 4, "degrees"),
    "PBLH": (("time",), 4, "m"),
    "release_lon": (("time",), 4, "Degrees_east"),
    "release_lat": (("time",), 4, "Degrees_north"),
    "particle_locations_n": (("height", "lon", "time"), 7, ""),
    "particle_locations_e": (("height", "lat", "time"), 7, ""),
    "particle_locations_s": (("height", "lon", "time"), 7, ""),
    "particle_locations_w": (("height", "lat", "time"), 7, ""),
}
_EDGE_NAMES = {"n": "N", "e": "E", "s": "S", "w": "W"}


def write_legacy_fp_file(ds, path, release_duration_hours=1 / 60, rounding="legacy"):
    """Write a footprint file in the legacy (2019) layout, with the same netCDF encoding.

    Args:
        ds: output of ``process_day_legacy`` (one or more days concatenated along ``time``).
        path: output file (overwritten if it exists).
        release_duration_hours: written into the time variable's ``period`` attribute, as legacy.
        rounding: "legacy" rounds each variable to its legacy ``least_significant_digit`` (lossy:
            footprints to multiples of 2**-17, which biases column totals LOW by 0-1.4%, X7), so
            the files can be compared with the 2019 ones exactly. "none" stores the full float32
            values (X7 adopted: use this for new footprints). The layout is the same either way.
    """
    if rounding not in ("legacy", "none"):
        raise ValueError(f"rounding must be 'legacy' or 'none', not {rounding!r}")
    import netCDF4

    times = ds["time"].to_index().to_pydatetime()
    reference = times.min().replace(month=1, day=1, hour=0, minute=0, second=0, microsecond=0)
    seconds = [(t - reference).total_seconds() for t in times]

    with netCDF4.Dataset(path, "w") as nc:
        nc.createDimension("time", ds.sizes["time"])
        nc.createDimension("lon", ds.sizes["lon"])
        nc.createDimension("lat", ds.sizes["lat"])
        nc.createDimension("lev", 1)
        nc.createDimension("height", ds.sizes["height"])
        nc.max_level = ds.attrs["max_level"]

        time = nc.createVariable("time", "d", ("time",))
        time[:] = seconds
        time.long_name = "time"
        time.standard_name = "time"
        time.units = f"seconds since {reference}"
        time.label = "left"
        time.period = f"{release_duration_hours} hours"
        time.comment = "time stamp corresponds to the beginning of each averaging period"
        time.calendar = "gregorian"

        for name, units in (("lon", "Degrees_east"), ("lat", "Degrees_north")):
            variable = nc.createVariable(name, "f", (name,))
            variable[:] = ds[name].values
            variable.units = units
        lev = nc.createVariable("lev", "S1", ("lev",))
        lev[:] = np.array(["column"])

        def digits(legacy_digits):
            return legacy_digits if rounding == "legacy" else None

        height = nc.createVariable(
            "height", "f", ("height",), zlib=True, least_significant_digit=digits(4)
        )
        height[:] = ds["height"].values

        for name, (dims, legacy_digits, units) in LEGACY_FILE_VARIABLES.items():
            variable = nc.createVariable(
                name, "f", dims, zlib=True, least_significant_digit=digits(legacy_digits)
            )
            variable[:] = ds[name].transpose(*dims).values
            variable.units = units
            if name.startswith("particle_locations_"):
                side = _EDGE_NAMES[name[-1]]
                variable.long_name = f"Fraction of total particles leaving domain ({side} side)"


def legacy_particles_from_fate(particles, z_units):
    """Turn a NAME 8.x fate-file table (``read_name_fate_file``) into the legacy particle table.

    The legacy processing expects heights in metres above ground in a column ``height_m_agl``
    (as in the ACRG particle files). Our NAME input writes the fate heights in m agl
    (``RunSettings.fate_z_coord``); anything else is refused rather than converted.
    """
    if z_units != "m agl":
        raise ValueError(f"Fate heights are in {z_units}, expected m agl")
    return particles.rename(columns={"z": "height_m_agl"})
