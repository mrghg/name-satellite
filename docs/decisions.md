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
| D11 | 2026-11-07 | **Clarity over perfect code.** Simple readable code; the workflow is clearly described in docstrings and `docs/workflow.md`. | Code must be checkable by scientists. |
| D12 | 2026-11-07 | `data/` is read-only. Legacy or experimental code runs on **copies** in scratch, never through symlinks; no `rm -rf` on variable paths. | Some reference files may be the last copies. |
| D13 | 2026-11-07 | **Particle-fate accounting for boundary sensitivities.** Every released particle is assigned to exactly one fate, and all fates are stored as fractions of the particles released, so they sum to exactly 1: (1) left through a lateral face of the dispersion domain → `N/E/S/W` histograms (along-edge cell × height), with the along-edge position clipped to the nearest output-grid cell; (2) left through the top → `top` histogram (lat × lon); (3) still inside at the end of the run → **only the fraction** is stored (per release and per column), not its location: if it is large the footprint should not be used, and if small the histograms can be rescaled with little error (decided with Matt, 2026-11-07); (4) anything else (NAME error flags) → a count that must be 0, else fail. Mean age per bin is kept (enough time information). How to treat (2) and (3) is left to the inversion; a legacy-style lateral-only renormalised version may be derived for compatibility. Height coordinate of the histograms: see open question 8. Legacy 2019 and 2025 variants are reproducible behind a switch (X14). | The inversion needs histograms that sum to 1 (background = histograms × global-model curtain). The legacy edge assignment double-counts corners and drops particles (X14), and ignores top exits and particles still inside (X6). |
| D14 | 2026-11-07 | Boundary-histogram heights are binned in **pressure** (fixed pressure bins, finer near the surface; NAME fate output requested in `Pa`). The legacy m-agl bins are reproducible behind the legacy switch. | Global-model curtains are on hybrid sigma-pressure levels, which convert to pressure simply (p = a + b·ps); retrievals and release layers are already in pressure; terrain is handled naturally (empty bins below the surface). To check: NAME 8.7.1 supports a `Pa` Z-coordinate for particle-information output. |
| D15 | 2026-11-09 | LOTUS2 jobs use `--account=acrg`, submitted from the sci servers. NAME runs in scratch; outputs to keep go to `~/acrg/mattrigby` (GWS nearly full: keep only processed footprints, inputs, logs). | Matt's choice; storage limits in `docs/jasmin.md`. |

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
4. ~~LOTUS2 account/QoS, scratch, output location, GitHub~~ Mostly resolved (2026-11-09): submit from
   sci servers; QoS/partitions in `docs/jasmin.md`; outputs to keep go to `~/acrg/mattrigby`
   (GWS nearly full), working files to scratch; GitHub `mrghg/name-satellite`; Slurm account `acrg`.
5. **GOSAT-2 footprint size** (legacy uses 10.5 km for GOSAT-1).
6. Keep legacy **binning** (means over 0.234° × 0.352° cells) as an option?
8. ~~**Height reference of boundary histograms.**~~ **Resolved by D14 (pressure).** Background: The legacy particle heights are **m agl** (checked
   2026-11-07: top exits sit at 19,011 m over the Altiplano and 19,012 m over the Atlantic, i.e. the
   19 km agl domain top; the ACRG writer outputs the particle's internal NAME height unconverted). The
   histogram height bins (0–20 km, 1 km) carry no reference, and the inversion's curtains may treat
   them as asl. That is fine over ocean or lowland edges, wrong over high terrain. Options for the new
   files: request the fate output in `m asl` or `Pa` (stock NAME converts) and bin on that. A
   decision on terrain handling is needed.
7. Details of the legacy **concentration unit-conversion bug**. Measured effect for Jan 2013: conc-based
   `fp` is about 1.42× (1.07–2.20×) the MR-based `fp` (docs/reference-data.md).
