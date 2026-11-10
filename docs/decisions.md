# Decisions

Newest at the bottom. Each entry gives the date, the decision and the reason. Decisions are made
with Matt Rigby unless noted.

| # | Date | Decision | Why |
|---|---|---|---|
| D1 | 2026-11-06 | Support both **time-integrated** and **hourly (HiTRes)** footprints. | Both are used by different inversions. |
| D2 | 2026-11-06 | GOSAT first; design for arbitrary satellites through a canonical sounding table. | GOSAT has the legacy reference; other instruments are planned. |
| D3 | 2026-11-06 | Target the **latest NAME** (8.7.1 on the NAME GWS). | Mk12 met needs 8.7.1; the legacy 7.2/8.4/8.6 ACRG builds are stale. |
| D4 | 2026-11-06 | Output contract: the netCDF written by legacy processing. Both `fp` (used by the current satellite-inversion code) and `srr` (newest inversions) are valid targets. **To be superseded** by the new footprint file (keyed by sounding ID, D18) once the group agrees the draft `docs/footprint-file-spec.md` (open question 13). | The satellite inversion team still runs the legacy `fp` code. |
| D5 | 2026-11-06 | Upper-level handling stays in the inversion (it reads `max_level`). The pipeline allows **any release ceiling below the NAME top**, specified in **pressure**. | Some systems only run to 5 or 10 km; retrieval level pressures vary per sounding. |
| D6 | 2026-11-06 | Use **mixing-ratio** field output, not air concentration. | The legacy concentration-based runs had unit-conversion errors (details to be found). |
| D7 | 2026-11-06 | First supported L2 input: **UoL OCPR proxy XCH4, GOSAT-1 v10.0 and GOSAT-2 v2.0** (CEDA, open access). ESA CCI v7.x is used only for legacy validation. | Newest products, covering 2009–2025; same file layout as legacy. |
| D8 | 2026-11-06 | Workflow: **reproduce legacy, then fix bugs as logged deviations, then add features**. | Avoid introducing silent changes. |
| D9 | 2026-11-06 | Tooling: uv, ruff, gitleaks (pre-commit), pytest; Python 3.12. Reference data stay out of git, verified by sha256. | |
| D10 | 2026-11-06 | Small runs on the sci node; NAME runs and heavy processing on LOTUS2. | JASMIN usage rules. |
| D11 | 2026-11-07 | **Clarity over perfect code.** Simple readable code; the workflow is clearly described in docstrings and `docs/workflow.md`. | Code must be checkable by scientists. |
| D12 | 2026-11-07 | `data/` is read-only. Legacy or experimental code runs on **copies** in scratch, never through symlinks; no `rm -rf` on variable paths. | Some reference files may be the last copies. |
| D13 | 2026-11-07 | **Particle-fate accounting for boundary sensitivities.** Every released particle is assigned to exactly one fate, and all fates are stored as fractions of the particles released, so they sum to exactly 1: (1) left through a lateral face of the dispersion domain → `N/E/S/W` histograms (along-edge cell × height), with the along-edge position clipped to the nearest output-grid cell; (2) left through the top → `top` histogram (lat × lon); (3) still inside at the end of the run → **only the fraction** is stored (per release and per column), not its location: if it is large the footprint should not be used, and if small the histograms can be rescaled with little error (decided with Matt, 2026-11-07); (4) anything else (NAME error flags) → a count that must be 0, else fail. Mean age per bin is kept (enough time information). How to treat (2) and (3) is left to the inversion; a legacy-style lateral-only renormalised version may be derived for compatibility. Height coordinate of the histograms: pressure (D14). Needs each particle exactly once, so it waits for X20 (duplicate fate records; interim: a patched NAME build that writes each particle's index). Legacy 2019 and 2025 variants are reproducible behind a switch (X14). | The inversion needs histograms that sum to 1 (background = histograms × global-model curtain). The legacy edge assignment double-counts corners and drops particles (X14), and ignores top exits and particles still inside (X6). |
| D14 | 2026-11-07 | Boundary-histogram heights are binned in **pressure** (fixed pressure bins, finer near the surface; NAME fate output requested in `Pa`). The legacy m-agl bins are reproducible behind the legacy switch. | Global-model curtains are on hybrid sigma-pressure levels, which convert to pressure simply (p = a + b·ps); retrievals and release layers are already in pressure; terrain is handled naturally (empty bins below the surface). How (2026-11-10): the fate output keeps heights in m agl (needed to find particles leaving through the 19 km agl top) and adds NAME's met columns (`RunSettings.fate_met`), which include the pressure at each particle's final position; every record has it. |
| D15 | 2026-11-09 | LOTUS2 jobs use `--account=acrg`, submitted from the sci servers. NAME runs in scratch; outputs to keep go to `~/acrg/mattrigby` (GWS nearly full: keep only processed footprints, inputs, logs). | Matt's choice; storage limits in `docs/jasmin.md`. |
| D16 | 2026-11-09 | For the NAME 8.7.1 validation runs against the legacy 2013 footprints, the bottom release layer uses the **legacy NAME surface pressure** (recovered from the legacy input files with `legacy.legacy_surface_pressure_hpa`), not the retrieval's own (X9 stays open for production runs). | Like-for-like comparison: differences then come only from NAME version/met handling, the integration-window fix (X18) and Monte Carlo noise. |
| D17 | 2026-11-10 | **GOSAT-1 and GOSAT-2 are separate observation streams.** Both UoL products are used (GOSAT-1 v10.0, 2009-04 to 2025-12; GOSAT-2 v2.0, 2019-02 to 2025-12), but never merged: separate observation files (and NAME runs/footprint files) per satellite, each recording its product and version. | Different instruments and different bias corrections (GOSAT-1: flat 16.655 ppb; GOSAT-2: 18.044 + 0.26 × SZA ppb), so the inversion can give each its own bias/offset. Matt, 2026-11-10. |
| D18 | 2026-11-10 | **Soundings are identified by their sounding ID** (the L2 exposure ID; for a legacy bin, the exposure ID of its earliest sounding, with all its exposure IDs kept in `exposure_id`). Our observation files are indexed by `sounding_id`, with time as an ordinary variable (`observations.to_obs`); release plans carry `sounding_id`. The legacy time-indexed layout is kept only for comparisons with legacy files. | GOSAT-1 v10 and GOSAT-2 v2 times are whole minutes, so times repeat (X22); time cannot be a key. Matt, 2026-11-10. Footprint files and the inversion's observation reading will need to use the same key. |

## Deferred options (evaluate after baseline validation)
- **(a) Weighted release:** one column field per sounding using Source Groups, with particles
  allocated in proportion to weight. Caveats: it ties the run to one retrieval version and
  species; weights can be negative or near zero; per-level diagnostics are lost.
- **(b) NAME NetCDF output:** route `N` with the NetCDF executable. Caveat: dense output size.
  Note: PR #23's "IA2" is the NAME II *text* format, not NetCDF.

## Open questions
Numbers are kept when a question is resolved (see the list below), so references stay valid.

1. **Bottom-layer surface pressure (X9).** Legacy used a separate NAME "SurfacePressure" run
   (`use_name_pressure=True`); its values differ from GOSAT p0 by −15 to +37 hPa. New runs
   (`prepare.py`) use the retrieval's p0, the only option for new data. First check (GOSAT-1 v10,
   1 Jan 2013, 46 soundings): the model surface pressure (release-point met) is a median 10.2 hPa
   higher than the retrieval's (−7.3 to +33.4), so the legacy offsets look like genuine model vs
   retrieval differences; no bottom release starts below the model ground. Proposal: adopt the
   retrieval's p0 and record the model's alongside.
3. **ACRG NAME code modification.** Stock 8.7.1 fate output suffices for footprints, but not for
   D13: it has no particle index, and duplicate records (X20) cannot be removed reliably. Agreed
   (2026-11-10): build a patched copy that adds the per-source particle index and fate code to the
   fate line. Still to decide: where the executable lives (proposed `~/acrg/mattrigby/NAME/`).
5. **GOSAT-2 footprint size** (legacy uses 10.5 km for GOSAT-1). The GOSAT-2 v2.0 files also say
   "10.5km x 10.5km at nadir", possibly copied from GOSAT-1 (GOSAT-2's is thought to be ~9.7 km).
   Matt is asking Leicester (2026-11-10), together with whether full-precision sounding times
   exist (X22).
6. **Binning** (legacy: means over 0.234° × 0.352° cells, X23): keep as an option or default?
   Without it, a month of GOSAT-1 v10 needs about 1.9 times the legacy NAME releases (Jan 2013:
   21,845 vs 11,526); the legacy binned uncertainty does not shrink with the number averaged.
9. **Which retrieval levels NAME models** (`max_level`, legacy 17 of 20). The paper's method
   (Ganesan et al. 2017) uses the prior above the NAME domain, but levels at 80 and 50 hPa are
   modelled although their releases mostly leave through the 19 km top (~4.5% of the column
   weight). Options: fewer levels, a pressure ceiling (D5, X10), or proper top-exit accounting
   (D13).
10. **Particle numbers for the upper levels.** Legacy releases 3/s above level 12 (180 particles);
    those releases are the noisiest (stage 1: 10–14% scatter between runs).
11. **Domains for new runs.** Matt is re-evaluating them (2026-11-10); no full month is run until
    then.
12. **Legacy hourly footprints (X3):** which fields format did the legacy hourly (HiTRes) runs
    write? If NAME III format, the hourly converter shifted those footprints by half a cell.
13. **New footprint file format:** the group's review of `docs/footprint-file-spec.md` (its own
    open questions: daily vs monthly files, pressure bins, top-exit grid, legacy-named copies,
    per-level footprints, who updates the inversion code).

## Resolved questions
2. **Release-point met** (the per-release met-extraction run exceeded NAME's grid limits):
   resolved 2026-11-09 by taking the met from the main run (`release_met.py`, X21).
4. **LOTUS2 account/QoS, scratch, outputs, GitHub:** resolved 2026-11-09 (D15, `docs/jasmin.md`).
7. **Legacy concentration unit-conversion bug:** resolved 2026-11-07 (X15): the 2018 processing
   converted with the air density at the release point, giving footprints about 1.44× too large.
8. **Height reference of boundary histograms:** resolved by D14 (pressure). The legacy particle
   heights are m agl (top exits sit at 19,011-19,013 m); the legacy height bins carry no reference.
