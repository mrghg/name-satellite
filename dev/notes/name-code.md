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

### Stage 1 (30-day) results, 2026-11-09
In the 30-day runs (stage 1, run B, `Sync? = Yes`, 6 min sync) the problem is much wider: 406/425
sources have extra records, 12,348 in total (2.2%). By level: ~35 extra per source for levels
1–12, ~3 for 13–16, ~61 for 17. Classified by end position, B has ~10,500 more lateral exits and
~1,650 more top exits than legacy, with the same number reaching the end of the run; the excess is
1–2% on every calendar day. No two records are identical, and near-identical pairs have random
height differences, so these are not the same particle written twice in the same state. Exit
times relative to the 6 min sync cycle: legacy never writes a record at 280–360 s into the cycle;
B writes ~12,700 there, plus 1,500 at 80–100 s.

Code reading so far: fate records are written for marked particles in `CalcParticleResultsSub`
(Case.F90 ~3509, sync-time `ReqInfoL` output), which then removes them
(`particleArrayRemoveParticle`, which doesn't compact the array). Newly released particles are
also written through `ReqInfoTPPInfos` (Case.F90 ~3880), which explains the doubled records for
releases above the domain top, but not the 30-day excess. Not resolved; a question for the NAME
team, with a minimal case (one source, a few days, fate output on).

## Point met output (release_met.py), checked 2026-11-09
- Fields can be output on an unstructured horizontal grid built from a `Locations:` set
  (`Horizontal Grids:` with `Name, dX, dY, Set of Locations`), and on a vertical grid given as an
  array (`Vertical Grids:` with `Z-Array`, values in an `Array:` block).
- A time grid given as a `T-Array` is read but can't be used for field output: `TInTGrid` etc.
  (GridAndDomain.F90 ~2640) assume T0 + i·dT, and the run aborts with "UNEXPECTED FATAL ERROR in
  MultiplyShortTime". Use one regular grid (nT = 1) per time instead.
- Limits: MaxArrays = 12, MaxArrayLength = 86, MaxLocationses = 5, MaxTGrids = 200; output
  groups default to 100 (`Max # Field Output Groups` in Main Options).
- Output file per output group `<group>_C1.txt`: fields-file layout with 15 header rows
  (request, quantity, units, ..., "Z = <value> <units>"), preliminary columns
  `T Index, Location Index, T, Location`, one row per location, 7 significant figures.
- A pressure (Pa) point below the model ground is reported at the ground: its "Pressure (Pa)" is
  exactly the surface pressure, and temperature etc. are the surface values.

## Output line length: at most ~2,218 field columns per output group (2026-11-10)
`GlobalParameters.F90`: `MaxOutputLineLength = 60000` (non-Compaq builds),
`OutputPrelimColumnWidth = 24`, `OutputFieldColumnWidth = 26`. `Output.F90` ~14970 (and ~15100 for
screen output) stops the run with "FATAL ERROR: too long an output line required for output group"
when `nPrelimCols*(PCW+1) + nFieldCols*(FCW+1) > MaxOutputLineLength`: 4 × 25 + n × 27 ≤ 60,000 gives
n ≤ 2,218 field columns per output group. The check runs only when the group is FIRST WRITTEN, which
for time-integrated output is the end of the run: stage 2 batch job 66323622 (2,720 sources in one
group) ran 8 h 21 min and then failed, with only the fate file written. `name_input` now splits the
field requests into groups of ≤ 2,000 (`Fields_grid1`, `Fields_grid1_2`, ...). Lesson: test new
run layouts with a 1-hour max-age run on the debug queue first.
