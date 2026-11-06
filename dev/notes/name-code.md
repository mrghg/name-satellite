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
