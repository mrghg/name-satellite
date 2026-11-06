# Reference data

Legacy files used to validate this pipeline. They are **not in git**. Put them under `data/`, or
point `NAME_SATELLITE_REFDATA` at a directory with the same layout. Paths and sha256 checksums are
registered in `src/name_satellite/refdata.py`, and `tests/test_refdata.py` checks them.

## GOSAT-BRAZIL, SOUTHAMERICA domain, 2013-01-01

There are **two independent NAME 7.2 runs** of this day with identical sources (25 soundings × 17
levels = 425 sources) and an identical observation file:

| | `gosat_brazil_201301/…_column_conc` (all of Jan 2013) | `gosat_brazil_20130101_mr_2019/…_column` (1 Jan only) |
|---|---|---|
| Run date | 17 June 2018 | 26 April 2019 |
| Origin | Matt's legacy zip archive | `NAME_processing/tests/files/LPDM/raw_output/Satellite_ByDay_MF` (commit 309e32e) |
| Field quantity | Air concentration (g s/m³) + BL-averaged fields | **Mixing ratio** |
| Dispersion domain | −93..−22 °E, −62.5..25 °N | output grid ±5° |
| Threads / mol. weight | 1 / 0 | 4 / 1.0 |
| Status | superseded (concentration unit problem) | **the run production `fp` came from** |

Shared by both: output grid 190 × 357 cells at 0.352° × 0.234°, origin (−91.329, −60.981);
domain Z max 19 000 m agl; 30-day back-runs; Mk6 global met; 60 s releases at 30/s (levels 1–12)
or 3/s (levels 13–17). Each particle file has 562,925 records: every particle is written exactly
once.

Notes:
- `Particle_files/Reprocessed` = `Original` with whitespace normalised only; the top-level
  `Particle_files` copy is the reprocessed one; `Original2` and `Particle_files_org` duplicate
  `Original`.
- `Met_Old` labels by level only (`_13_`); `Met` uses point+level IDs (`_02107_`).
- Observations: ESA CCI OCPR **v7.2**, quality-filtered, NAME-surface-pressure "level_order"
  filtered, binned to 0.234° × 0.352°. v7.2 is not on CEDA, so this file *is* the ground truth.

### January 2013, concentration run (`gosat_brazil_201301/`)
Downloaded from Matt's archive on 2026-11-06. It supersedes the earlier 1-day copy: all 864 files of
that copy were byte-identical. Three stray Met folders (2010-05-21, 2014-06-08, 2015-03-31, caught by
the `*201301*` wildcard inside coordinates) were removed. Checks made (`tests/test_refdata.py`):
- all 31 runs logged "Run completed successfully" and have empty error files;
- per day: sources = soundings × 17; particle records = Σ(1801 or 181 per source) exactly; met files
  and field columns = number of sources;
- 678 soundings in total; observation times identical to the production file;
- sha256 of all 435 non-met files in `dev/manifests/gosat_brazil_201301_conc.sha256`; the met folders
  are checked by count. `Met_Old` is incomplete for 2013-01-07 (323 of 1088; its level-only names
  collide).

`Processed_Fields_files/GOSAT-BRAZIL-column_SOUTHAMERICA_201301.nc` is the **concentration-run**
`fp` (internal `created` 2018-09-14). Its release met is identical to production, but its footprints
are much larger: total-`fp` ratio conc/MR has median **1.42**, 5–95% 1.07–2.20, range 0.98–3.0. That
is far beyond Monte Carlo noise, and is presumably the legacy concentration unit-conversion error
(decision D6).

## Production processed file

`GOSAT-BRAZIL-column_SOUTHAMERICA_201301.nc`: January 2013, 31 days, 678 soundings, `fp` in
(mol/mol)/(mol/m²/s), `max_level=17`, created 2019-05-04 by rt17603. Its 2013-01-01 slice is
**bit-identical** to the NAME_processing benchmark (`fp`, met, release positions and
histograms), so it was made from the 2019 mixing-ratio run.

## Validation chain

1. Obs file (+ NAME surface pressure) → NAME source definitions: reproduced exactly
   (`dev/scripts/reproduce_legacy_sources.py`).
2. 2019 MR raw output → benchmark/production `fp`: target for the processing port (exact).
3. 2018 conc run vs 2019 MR run: Monte Carlo noise (1 Jan), and the unit-conversion error (whole
   month, via the two processed files).

## TODO
- ☑ January 2013 raw output of the **2018 concentration run** (added 2026-11-06).
- ☐ January 2013 raw output of the **2019 mixing-ratio run** (only 1 Jan so far), for month-scale
  exact checks against the production file. Maybe on the other netCDF/archive server.
- ☐ A recent (Mk11/Mk12-era) successful legacy run, if one exists.

## Other known locations (JASMIN)
- Processed `fp` examples (older format, 2019): `/gws/ssde/j25b/acrg/elenafi/example_footprints/`
  (e.g. GOSAT-BRAZIL 201801, CHINA 201605, INDIA 201503, SAHARA 2014–16).
- NAME 8.4 inputs from a failed Oct 2024 run (Alaska, 2018-07-01, 1853 sources):
  `/gws/ssde/j25b/acrg/adanjou/NAME/Results/ALASKA_GOSAT-ALASKA_column`.
- L2 data: see `dev/notes/gosat-l2-and-preprocessing.md`.
