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

## Current products and the new reader (checked 2026-11-09)
Files used: GOSAT-1 v10.0 for 2013-01-01/02/07, GOSAT-2 v2.0 for 2020-01-01/02 (dap.ceda.ac.uk,
~1 MB each), ESA CCI v7.0 for 2013-01-01/07 (/neodc/esacci), all copied to
/work/scratch-pw5/mattrigby/name-satellite/l2_20261109.
- **Times are whole minutes in v10 and GOSAT-2 v2** (X22): v10 time = v7.0 time floored to the
  minute for all 2,024 common soundings on 2013-01-01. 2,179 soundings share 305 distinct times.
  The exposure ID's embedded time is not the measurement time either (e.g. ID ...14:58:02 for a
  sounding at 14:58:50).
- The legacy `gosat_process_file` crashes on v10/GOSAT-2 (binning: "cannot reindex ... duplicate
  values"); it never processed these products. With today's xarray its own file writer also fails
  (exposure IDs as bytes), so `dev/scripts/run_legacy_gosat.py` saves the dataset before writing.
- Running the legacy code: jaspy 3.12 lacks `rich`, `iso3166`, `openghg_defs`; they were installed
  with `uv pip install --target` into /work/scratch-pw5/mattrigby/name-satellite/legacy_pydeps_20261109.
- GOSAT-2 v2.0 declares `_FillValue = -9999.99` on pressure levels, so xarray reads them as NaN and
  the legacy bad-pressure filter matches nothing (X24); 6 good-quality soundings on 2020-01-01 have an
  undefined surface level. ESA CCI v7.0 has NaN levels only in quality-flagged soundings.
- Some soundings over high ground (Tibet, Antarctica) have level 0 (retrieved surface pressure)
  LOWER than level 1; 2 such good-quality soundings on 2013-01-07 (v7.0). `level_order_mask` drops them.
- v10 attributes note a bias correction: "GOSAT-1 flat offset 16.655 ppb subtracted ...".
  `ch4_profile_apriori` units are "1e-9" though its description says ppm (values are ppb).
- Counts, 2013-01-01 v10: 2,179 soundings, 863 good; binned in the GOSAT-BRAZIL box: 23 (v7.0: 25;
  the 2013 production v7.2 file: 25).
- Our `observations.py` matches the legacy preprocessing exactly on v7.0 (1 and 7 Jan: 50 / 92
  unbinned, 25 / 64 binned; all values identical; dlat/dlon to 1e-7 because legacy stored float32).
