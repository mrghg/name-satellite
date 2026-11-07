# Deviations from the legacy workflow

Every change that makes our numbers differ from the legacy pipeline must be listed here, with
evidence and a test that quantifies its effect. Until an entry is marked **Adopted**, the code
reproduces legacy behaviour (possibly behind a `legacy=True` switch).

Status: **Candidate** (suspected, not yet confirmed) → **Confirmed** (demonstrated with a test)
→ **Adopted** (fixed in default behaviour) / **Rejected**.

| ID | Status | Area | Legacy behaviour | Proposed change | Evidence / test |
|---|---|---|---|---|---|
| X1 | Candidate | Met version | `get_Mk_Global` in `BackRuns_satellite.scr` stops at Mk10, so all post-2017 dates use Mk10 (whose archive ends 2022-05-04 12Z). PR #29 adds Mk11 with a boundary that differs from the site script (2022-05-05 vs 2022-05-03). | Choose the met version from actual archive coverage. | `dev/notes/legacy-run-scripts.md` |
| X2 | Candidate | Run length | MaxAge reduced by hand around the Mk6→7 and Mk8→9 changeovers in `Satellite_launch.scr`. | Derive it from met coverage; record it in the output. | same |
| X3 | Candidate | Processing grid | The hourly converter adds +0.5 cell to the grid origin, while the text path treats the origin as the cell centre (possible half-cell shift). | Check against NAME's grid definition; fix. | `dev/notes/legacy-processing.md` |
| X4 | Candidate | Matching | Fields, met, particles and obs are aligned by order only. | Match by source/sounding ID. | same |
| X5 | Candidate | Particle IDs | Particle Id → (time, level) uses a running counter, assuming contiguous IDs. | Map via source index/name. | same |
| X6 | Candidate | Histograms | Boundary histograms are normalised by edge particles only; top exits are excluded. | Decide and document. | same |
| X7 | Candidate | Precision | `srr` is stored with `least_significant_digit=5` (lossy). | Lossless or justified precision. | same |
| X8 | Candidate | Metadata | `release_height` is in Pa but labelled m. | Fix the units attribute. | same |
| X9 | Candidate | Surface pressure | Bottom release layer is anchored to the surface pressure from a separate NAME run (differs from GOSAT p0 by −15 to +37 hPa for 2013-01-01). | Use the retrieval p0 and log the UM comparison (decision pending). | `dev/scripts/reproduce_legacy_sources.py` |
| X10 | Candidate | Upper levels | Level 17 (about 65–30 hPa) is released partly above the 19 km domain top; 89% of its particles die at age 0. | Pressure ceiling checked against the domain top (D5). | `dev/notes/legacy-processing.md` |
| X11 | Candidate | Error handling | Many silent `return None`s in processing. | Raise. | — |
| X12 | Candidate | Binning time | `ds_check_internal_unique` adds *random* µs offsets to duplicate times (not reproducible). | Deterministic IDs instead. | `dev/notes/gosat-l2-and-preprocessing.md` |
| X13 | Candidate | Layer edge case | `midpoint_bounds` top-layer fix-up assigns a full-length array to a subset (would fail if only some soundings hit the 0 hPa limit). | Vectorised fix. | same |
| X14 | Candidate | Histograms (versions) | 2019 production histograms differ slightly from today's legacy code (same totals per sounding; per-edge/height redistribution). Changed between May 2019 and Sept 2021. | Decide which version to reproduce; find the commit. | docs/reference-data.md (legacy processing check) |
| X15 | Candidate | g s/m³ conversion | Conc fields are converted with a fixed surface P/T = 345 (`use_surface_conditions=True`) instead of met at the release. The 2018 processing code gave footprints about 1.44× too large. | Use mixing-ratio output only (D6); if concentration is ever needed, use the release-point met. | same |
| X16 | Candidate | Histogram gaps | When no particles reach an edge, processing logs "Copying lower level/previous time step" and fills the histogram from a neighbour. | Check what is copied and whether it should be. | `run_legacy_process` logs |
