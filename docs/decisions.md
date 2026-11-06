# Decisions

Newest at the bottom. Each entry gives the date, the decision and the reason. Decisions are made
with Matt Rigby unless noted.

| # | Date | Decision | Why |
|---|---|---|---|
| D1 | 2026-11-06 | Support both **time-integrated** and **hourly (HiTRes)** footprints. | Both are used by different inversions. |
| D2 | 2026-11-06 | GOSAT first; design for arbitrary satellites through a canonical sounding table. | GOSAT has the legacy reference; other instruments are planned. |
| D3 | 2026-11-06 | Target the **latest NAME** (8.7.1 on the NAME GWS). | Mk12 met needs 8.7.1; the legacy 7.2/8.4/8.6 ACRG builds are stale. |
| D4 | 2026-11-06 | Output contract: the netCDF written by legacy processing. Both `fp` (used by the current satellite-inversion code) and `srr` (newest inversions) are valid targets. | The satellite inversion team still runs the legacy `fp` code. |
| D5 | 2026-11-06 | Upper-level handling stays in the inversion (it reads `max_level`). The pipeline allows **any release ceiling below the NAME top**, specified in **pressure**. | Some systems only run to 5 or 10 km; retrieval level pressures vary per sounding. |
| D6 | 2026-11-06 | Use **mixing-ratio** field output, not air concentration. | The legacy concentration-based runs had unit-conversion errors (details to be found). |
| D7 | 2026-11-06 | First supported L2 input: **UoL OCPR proxy XCH4, GOSAT-1 v10.0 and GOSAT-2 v2.0** (CEDA, open access). ESA CCI v7.x is used only for legacy validation. | Newest products, covering 2009–2025; same file layout as legacy. |
| D8 | 2026-11-06 | Workflow: **reproduce legacy, then fix bugs as logged deviations, then add features**. | Avoid introducing silent changes. |
| D9 | 2026-11-06 | Tooling: uv, ruff, gitleaks (pre-commit), pytest; Python 3.12. Reference data stay out of git, verified by sha256. | |
| D10 | 2026-11-06 | Small runs on the sci node; NAME runs and heavy processing on LOTUS2. | JASMIN usage rules. |

## Deferred options (evaluate after baseline validation)
- **(a) Weighted release:** one column field per sounding using Source Groups, with particles
  allocated in proportion to weight. Caveats: it ties the run to one retrieval version and
  species; weights can be negative or near zero; per-level diagnostics are lost.
- **(b) NAME NetCDF output:** route `N` with the NetCDF executable. Caveat: dense output size.
  Note: PR #23's "IA2" is the NAME II *text* format, not NetCDF.

## Open questions
1. **Bottom-layer surface pressure.** Legacy used a separate NAME "SurfacePressure" run
   (`use_name_pressure=True`); its inferred values differ from GOSAT p0 by −15 to +37 hPa (see
   `dev/notes/gosat-l2-and-preprocessing.md`). Proposal: use the retrieval's p0, compare it with
   the UM surface pressure and log the difference.
2. **Release-point met.** The per-release met-extraction run exceeds NAME's grid limits. Options:
   per-particle met output; reading UM/NAME met directly.
3. Is any **ACRG NAME code modification** still needed? Stock 8.7.1 fate output appears sufficient.
4. LOTUS2 **account/QoS**, scratch location, output location on the GWS, and the GitHub repo
   location.
5. **GOSAT-2 footprint size** (legacy uses 10.5 km for GOSAT-1).
6. Keep legacy **binning** (means over 0.234° × 0.352° cells) as an option?
7. Details of the legacy **concentration unit-conversion bug**. Measured effect for Jan 2013: conc-based
   `fp` is about 1.42× (1.07–2.20×) the MR-based `fp` (docs/reference-data.md).
