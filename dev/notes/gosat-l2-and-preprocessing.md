# GOSAT L2 data and legacy preprocessing (acrg.satellite)

Checked 2026-11-06.

## Products
| Product | Where | Coverage | Access |
|---|---|---|---|
| UoL GOSAT-1 Proxy XCH4 **v10.0** | `https://dap.ceda.ac.uk/neodc/gosat/data/ch4/UoL_gosat1_v10/CH4_GOSAT_OCPR/v10.0/YYYY/UoL-GHG-L2-CH4-GOSAT-OCPR-YYYYMMDD-fv10.0.nc` | 2009-04-23 → 2025-12-18 | Open (OGL v3); DOI 10.5285/94d6572a308e477c86e0ce7c4d771311 |
| UoL GOSAT-2 Proxy XCH4 **v2.0** | `…/UoL_gosat2_v2/CH4_GOSAT2_OCPR/v2.0/YYYY/UoL-GHG-L2-CH4-GOSAT2-OCPR-YYYYMMDD-fv2.0.nc` | 2019-02-05 → 2025-12-31 | Open (OGL v3); DOI 10.5285/8209b9acb92e4ba69188c9c8ff7a3b76 |
| ESA CCI OCPR v7.0 | `/neodc/esacci/ghg/data/crdp_4/GOSAT/CH4_GOS_OCPR/v7.0/` (mounted) | 2009-04 → 2015 | Open |
| ESA CCI OCPR v7.2 | not on CEDA | — | used by the legacy 2013 run |

- The JASMIN mount `/neodc/gosat` gives *permission denied* (group `byacl`) even though the data
  are open. Download over HTTPS from `dap.ceda.ac.uk` (works without login), or ask the JASMIN
  helpdesk.
- Layout (v7.0, v10.0, GOSAT-2 v2.0 all the same): dims `n` (soundings) × `m` (20 levels);
  variables used: `xch4, xch4_uncertainty, latitude, longitude, time, pressure_levels (hPa),
  pressure_weight, xch4_averaging_kernel, ch4_profile_apriori, exposure_id (char 22),
  xch4_quality_flag (0 = good), retr_flag (0 = land, 1 = glint)`. v10 and GOSAT-2 add e.g.
  `xch4_no_correction`.
- Pressure levels: the top five are fixed (80, 50, 10, 1, 0.1 hPa). Below that they vary per
  sounding (e.g. GOSAT-2: 945, 917, 864 … 177 hPa). Level 0 is the surface. Pressure weights
  sum to 1, with about half weight on level 0 (level-based quadrature).
- 2013-01-01: v7.0 has 2109 retrievals / 791 good; v10.0 has 2179 / 863. **Not the same data
  as v7.2.**

## Legacy preprocessing: `~/code/acrg/acrg/satellite/{gosat,common}.py`
`gosat_process(...)` → per file `gosat_process_file`:
1. `add_coords` (n, m → time, lev), sort by time.
2. Filters: quality flag (`xch4_quality_flag == 0`), bad pressure (any `-9999.99`), optional mode,
   optional NAME-surface-pressure filters (`level_order`, `cutoff`, `dpressure_range`).
3. Lat/lon bounds (or domain); optional binning (`create_zbin` + `mean_ds`: means within lat/lon bins
   over the whole file, time = midpoint, uncertainties combined, `exposure_id` concatenated with
   commas, `dlat/dlon` = footprint + spread). `ds_check_internal_unique` then adds random µs to
   duplicate times (X12).
4. `output()` writes the obs netCDF `gosat-fts_gosat_YYYYMMDD[-NNN]_ch4-column.nc` (dims time, lev, id).
5. `output_name()` writes the NAME CSV `SITE_YYYYMMDD[-A].csv`:
   - columns `ID_Level, Time, x, dx, y, dy, z, dz`;
   - ID = `PPP` + `LL` when `file_per_day`;
   - z/dz in **Pa** from `define_pressure_levels` (×100 from hPa);
   - levels `[:max_name_level=17]`;
   - dx/dy = haversine size of a **10.5 km** footprint (full width), e.g. dlat 0.0944°.

The NAME CSV and obs file are paired **only by row order** within a date.

### Release layers: `common.midpoint_bounds(set_as_edges=False)`
With p = level pressures and Δᵢ = pᵢ − pᵢ₊₁:
- level 0: dp₀ = Δ₀/2, centre = p₀ − dp₀/2. The layer runs from the surface to the midpoint of p₀ and p₁.
- levels 1..n−2: dpᵢ = (Δᵢ₋₁ + Δᵢ)/2, centre = pᵢ + Δᵢ₋₁/2 − dpᵢ/2. Boundaries are the midpoints.
- top: dp = Δ_last, with a fix-up if the layer would go below 0 hPa (edge-case bug X13).
- With `use_name_pressure=True`, p₀ is replaced by the NAME surface pressure before this.

**Verified:** all 425 sources in the 2013-01-01 legacy input are reproduced from the obs file to
float precision (Z rel. diff 7e-8, dZ 4e-7, X/Y 5e-11) **only** with the NAME surface pressure for
p₀. The inferred NAME p₀ differs from the GOSAT p₀ by −15 to +37 hPa (mostly higher; the first
soundings are ≈1013–1017 hPa). That is suspicious (sea-level pressure? a stale nearest-day template?). Script:
`dev/scripts/reproduce_legacy_sources.py`. The surface-pressure files came from a separate NAME
run (`Pressure_C1_*.txt`, `$DATA_PATH/LPDM/surface_pressure/<domain>`, nearest day within 31 days).
