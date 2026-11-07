# Workflow

This page describes the end-to-end workflow in plain terms. It describes the **target** design;
steps not yet built are marked *(planned)*. Each step names the module that does it, so you can go
from here to the code. Each module's docstring describes the same step in more detail.

```
 Satellite L2 files ─► 1. Read soundings ─► 2. Plan releases ─► 3. Write NAME inputs
                                                                        │
         footprint netCDF ◄─ 6. Combine into column ◄─ 5. Read NAME output ◄─ 4. Run NAME (LOTUS2)
```

## 1. Read soundings *(planned)*
**In:** one day of satellite Level-2 files (first: University of Leicester GOSAT / GOSAT-2 CH4
proxy). **Out:** a table with one row per sounding: ID, time, latitude/longitude, footprint size,
and per-level pressure (hPa), pressure weight, averaging kernel and a priori.
Quality filters (and optional legacy-style binning onto the NAME grid) are applied here. The same
table is written as the observation file that the inversion reads.

## 2. Plan releases *(in progress: `name_satellite/releases.py`)*
**In:** the sounding table and a pressure ceiling. **Out:** a release plan: for each sounding and
each retrieval level below the ceiling, a vertical layer (pressure top and bottom, Pa), a
horizontal box (the footprint size), a 60 s release time and a particle rate.
Layers run between the midpoints of neighbouring retrieval levels. The lowest layer starts at the
surface. Legacy details: `dev/notes/gosat-l2-and-preprocessing.md`.

Done: `release_layers_pa` (layer centre and thickness, Pa) and `particles_per_second` (legacy
rate rule). They reproduce the vertical layers, positions, times and rates of all 11,526 legacy
sources for January 2013 (`tests/test_releases.py`). Still to do: pressure ceiling, footprint size,
and assembling a release plan table.

## 3. Write NAME inputs *(planned)*
**In:** the release plan, domain settings and the met version available for the dates. **Out:** a
NAME input file per day: one source and one output request per release, plus particle-fate output.

## 4. Run NAME *(planned)*
Each day is a SLURM job on LOTUS2. NAME runs backwards in time for up to 30 days. It writes, for
each release, a footprint field (time-integrated or hourly) and a record of where and when each
particle left the domain.

## 5. Read NAME output *(in progress: `name_satellite/name_output.py`, `legacy.py`)*
**In:** NAME fields and particle files. **Out:** for each release, a footprint (sensitivity of the
mole fraction at the release to surface fluxes, (mol/mol)/(mol/m²/s)), counts of particles leaving
through each boundary (for boundary conditions) and the met at the release point.
Everything is matched by release ID, never by order.

Done: `read_fields` (NAME text fields → array over source × lat × lon, keyed by source name;
missing all-zero rows restored), `read_acrg_particle_file` (legacy particle fates), and legacy
release-point met (`legacy.read_legacy_met_folder`). Tested on every day of January 2013. Still to
do: hourly fields, and the stock NAME 8.7.1 fate-output format.

Particle fates (decision D13): every released particle ends in exactly one category. It either
left through the north, east, south or west face, left through the top, or was still inside at the
end of the run. Edge and top exits are stored as histograms of *fractions of the particles released*;
for particles still inside only their fraction is stored. All categories together sum to exactly 1. The inversion multiplies these by a global model's mole
fractions to get the background, and decides how to treat top exits and particles still inside.

## 6. Combine into column footprints *(planned)*
**In:** the per-level footprints and the sounding table. **Out:** one footprint per sounding,
Σ over levels of (pressure weight × averaging kernel × level footprint), written in the netCDF
format the inversion reads (legacy `fp` format and/or `srr`), with `max_level` recorded.

## Validation
Each step is first checked against the legacy pipeline's output for GOSAT-BRAZIL, January 2013
(see [reference-data.md](reference-data.md)), before any deliberate improvements
([deviations.md](deviations.md)).
