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
| X6 | Adopted (D13, to implement) | Histograms | Boundary histograms are normalised by edge particles only; top exits and particles still inside are excluded. | Account for every particle (D13). | same |
| X7 | Candidate | Precision | `srr` is stored with `least_significant_digit=5` (lossy). | Lossless or justified precision. | same |
| X8 | Candidate | Metadata | `release_height` is in Pa but labelled m. | Fix the units attribute. | same |
| X9 | Confirmed | Surface pressure | Bottom release layer is anchored to the surface pressure from a separate NAME run (differs from GOSAT p0 by −15 to +37 hPa for 2013-01-01). **Confirmed consequence:** in January 2013, 22 of 11,526 releases (21 of 678 soundings) are centred below the model's ground: NAME's met at the release reports a pressure 0.1–23 hPa lower than requested (levels 1 and 2 only), so part of the release layer was underground. | Use the retrieval p0 and log the UM comparison (decision pending). | `dev/scripts/reproduce_legacy_sources.py` |
| X10 | Candidate | Upper levels | Level 17 (about 65–30 hPa) is released partly above the 19 km domain top; 89% of its particles die at age 0. | Pressure ceiling checked against the domain top (D5). | `dev/notes/legacy-processing.md` |
| X11 | Candidate | Error handling | Many silent `return None`s in processing. | Raise. | — |
| X12 | Candidate | Binning time | `ds_check_internal_unique` adds *random* µs offsets to duplicate times (not reproducible). | Deterministic IDs instead. | `dev/notes/gosat-l2-and-preprocessing.md` |
| X13 | Candidate | Layer edge case | `midpoint_bounds` top-layer fix-up assigns a full-length array to a subset (would fail if only some soundings hit the 0 hPa limit). | Vectorised fix. | same |
| X14 | Adopted (D13, to implement) | Boundary histograms (edge assignment) | Particles are put on an edge if they lie beyond the output grid's outer half-cell, then binned along the edge with `np.histogram2d` over the grid cells. Because the 2019 MR run's dispersion domain is 5° larger than the grid, exiting particles sit up to 5° outside it: (i) 8,011 corner particles (1.5%) count on two edges; (ii) particles beyond the *lower* end of an edge are silently dropped in all versions (12,645); (iii) beyond the *upper* end they were dropped until NAME_processing 094d8e6 (2025-05-15) and are now put in the last cell (2,765). Histograms are then renormalised over what remains. The 2019 production files predate (iii); today's legacy code includes the asymmetric fix. | Assign each particle to exactly one exit boundary (the dispersion-domain face it crossed, or the top), clip its along-edge position to the nearest grid cell, normalise by particles released. Reproduce both legacy variants behind a switch for comparison. | dev/notes/legacy-processing.md (X14) |
| X15 | Confirmed | g s/m³ conversion | Converting air concentration (g s/m³) to mixing ratio needs the air density where the field is measured (0–40 m). The 2018 code (acrg d8a5e33d `convert_to_ppt`) used P/(RT) at the **release point**, which is much smaller aloft, so upper-level footprints were inflated (Jan 2013: about 1.44× too large overall, up to 3×). 2019 (75793919) kept this for g s/m³ but used the correct ppm-s path for mixing-ratio runs; later code switched g s/m³ to a fixed surface P/T = 345. | Use mixing-ratio output only (D6); if concentration is ever needed, use the near-surface air density. | docs/reference-data.md; acrg history |
| X16 | Candidate | Histogram gaps | When no particles reach an edge, processing logs "Copying lower level/previous time step" and fills the histogram from a neighbour. | Check what is copied and whether it should be. | `run_legacy_process` logs |
| X17 | Confirmed | Met units | The 2019 files store the release-point temperature in **°C** but label it `K` (`write_netcdf` units vs `read_met` °C). Today's legacy code converts to K. | Store K, labelled K. | 1 Jan 2013 golden test (production = met °C) |
