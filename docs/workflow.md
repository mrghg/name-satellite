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

## 2. Plan releases *(planned)*
**In:** the sounding table and a pressure ceiling. **Out:** a release plan: for each sounding and
each retrieval level below the ceiling, a vertical layer (pressure top and bottom, Pa), a
horizontal box (the footprint size), a 60 s release time and a particle rate.
Layers run between the midpoints of neighbouring retrieval levels. The lowest layer starts at the
surface. Legacy details: `dev/notes/gosat-l2-and-preprocessing.md`.

## 3. Write NAME inputs *(planned)*
**In:** the release plan, domain settings and the met version available for the dates. **Out:** a
NAME input file per day: one source and one output request per release, plus particle-fate output.

## 4. Run NAME *(planned)*
Each day is a SLURM job on LOTUS2. NAME runs backwards in time for up to 30 days. It writes, for
each release, a footprint field (time-integrated or hourly) and a record of where and when each
particle left the domain.

## 5. Read NAME output *(planned)*
**In:** NAME fields and particle files. **Out:** for each release, a footprint (sensitivity of the
mole fraction at the release to surface fluxes, (mol/mol)/(mol/m²/s)), counts of particles leaving
through each boundary (for boundary conditions) and the met at the release point.
Everything is matched by release ID, never by order.

## 6. Combine into column footprints *(planned)*
**In:** the per-level footprints and the sounding table. **Out:** one footprint per sounding,
Σ over levels of (pressure weight × averaging kernel × level footprint), written in the netCDF
format the inversion reads (legacy `fp` format and/or `srr`), with `max_level` recorded.

## Validation
Each step is first checked against the legacy pipeline's output for GOSAT-BRAZIL, January 2013
(see [reference-data.md](reference-data.md)), before any deliberate improvements
([deviations.md](deviations.md)).
