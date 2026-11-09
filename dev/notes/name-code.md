# NAME source code facts

From `~/name/code/NAME-8.7.1/Code_NameIII` (stock) and the ACRG build
`/gws/ssde/j25b/acrg/NAME/NAMEIII_v8_6_withUKV`. Checked 2026-11-06.

## Particle fate output (stock)
- Particles are written to Fate Info output only when `Marked ≠ 0` (Case.F90:6119).
- Marking codes (Particle.F90:103-109): M_Error=1, M_NotReleased=2, M_Transformed=4,
  **M_Location=8** (left the domain, including through the top), **M_Time=16** (max travel time / end of run).
- Columns (Output.F90 ~17230): `Source, Travel Time (hr), X, Y, Z`, comma-separated. Source is the
  integer source index.
- In the 2013 legacy runs every released particle appears exactly once.

## ACRG "AJM" modifications (v8.6 withUKV vs stock 8.6)
- Particle output reformatted to whitespace `X(F7.2) Y(F7.2) Z(F6.0) Age(F7.2) iSource`, with no
  PPInfo header. Same information as stock, at lower precision.
- `Dt = Min(Dt, Flow%MaxDT)` (UKV only).
- MaxHGrids/ZGrids/TGrids 70/20/200 → 500; MaxNWPMetDefns etc. → 32; OutputFieldColumnWidth 26 → 30.
- Conclusion: stock fate output is sufficient (needs a new parser). The raised grid limits mattered
  only for the per-release met-extraction run, which doesn't scale anyway.

## Limits (GlobalParameters.F90, stock 8.7.1)
MaxHGrids=70, MaxZGrids=20, MaxTGrids=200, MaxPPInfoReqs=10, OutputFieldColumnWidth=26,
MaxFieldReqsPerReqInfo=3000.

## Output routes and formats
- `FieldReq%NetCDF = Scan(OutputRoute,'Nn')` (Output.F90:4044); needs the NetCDF executable
  (check ~4127); incompatible with NAME II format `2` (4120).
- Output Format letters: I = indices, A = align, Z = zero lines suppressed, 2 = NAME II text.
- `NetCDF Output Options:` Use Shuffle?, Use Fletcher32?, Deflate Level (default 4).
- Source Groups (Input.F90:9303): a field request can be restricted by Source **or** Source Group
  (Output.F90:3646-3661). This enables the weighted-release option.

## Executables
8.7.1: `nameiii_64bit_par.exe`, `nameiii_eccodes_NetCDF_64bit_par.exe`. 8.6 also had
`nameiii_NetCDF_64bit_par.exe`. Mk12 met definitions exist only in 8.7.1.

## Threading benchmark (2026-11-09; 1 Jan 2013, 425 releases, 24 h, Mk6, debug queue)
| Threads | Total | Particle loop | Other | Speed-up | Efficiency |
|---|---|---|---|---|---|
| 1 | 813 s | 639 s | 174 s | 1× | 100% |
| 4 | 253 s | 162 s | 91 s | 3.2× | 80% |
| 8 | 179 s | 92 s | 87 s | 4.5× | 57% |

The particle loop scales almost linearly; the rest (met reading, output) barely changes. With
`Random Seed = Fixed`, fields and fate files are **bit-identical** for 1, 4 and 8 threads.
(NAME's timer summary is in `<input>Timing.txt`; `sacct` TotalCPU shows 0 under mpirun.)

## Duplicate fate records for particles released above the domain top (open; ask the NAME team)
In the benchmark, 22 of 425 sources have more fate records than particles released (60 × rate
+ 1), 1,616 in total. All the extras are particles released ABOVE the 19 km agl domain top and
recorded within 0.1 h. Level-17 releases of soundings 2–6 (release ends 15:01:27–15:02:13) have
exactly twice the particles (362 vs 181); soundings 7–10 have progressively fewer extras; sounding 1
(14:58:50) and most late soundings (16:34+) have none. The pattern depends on release time,
presumably relative to NAME's internal time steps. The legacy ACRG 7.2 files always had exactly
60 × rate + 1 records.

The fate file has no column giving the reason a particle stopped (left the domain vs reached the
end of the run vs not released), so we can't tell duplicates apart. This doesn't affect the
fields (these particles are outside the domain) or the legacy edge histograms (they're inside the
lateral bounds, so never counted), but it does matter for D13 (fractions per fate). Options: ask
the NAME team; keep releases below the domain top (D5 ceiling); count top exits per source and
cap at the number released.

`M_NotReleased` (Case.F90 ~8333) is for particles that reflect off the ground at release when the
source has `No Reflect? = Yes`, so it is not this effect.
