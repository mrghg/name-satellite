# NAME 8.7.1 input files: syntax and traps

Observed in our runs (Nov 2026). Working generator: `~/code/name-satellite/src/name_satellite/name_input.py`.

## General shape
- Plain text made of blocks: a heading line ending in `:`, a header row of column names, then one
  row per item, comma-separated. NAME reads columns **by heading**, so column order is free (rows
  follow their own header's order). Unknown columns give a warning unless the heading starts
  with `!` (B01).
- **A blank line ends a block**: never put blank lines inside one. Text outside blocks is a
  comment; start it with `!` or `#` (B01 Appendix).
- Every row must give one item per heading, so **a row whose last item is blank ends with a
  comma** (the separator); other trailing commas are optional (B01 Appendix).
- `Restart File Options:` and `Multiple Case Options:` blocks were needed in our 8.7.1 runs even
  when unused (our observation; the docs only say they are usually blank / 1). Copy them from the
  NAME team's 8.7.1 demo templates.
- `Main Options`: `Absolute or Relative Time?, Fixed Met?, Flat Earth?, Run Name, Random Seed,
  Max # Sources, Max # Field Reqs, Max # Field Output Groups, Backwards?`.
  - `Max # Sources` defaults to 100; raise it (no hard maximum documented).
  - `Max # Field Reqs` (default 650) must be well above the number of requests: in our runs
    mixing-ratio output needed extra internal fields (not documented). 20,000 works for ~3,000.
  - `Max # Field Output Groups` defaults to 100; raise it if you have many output groups.
  - `Same Results With/Without Update on Demand?` forces consistent sync steps (A01 §1.2.7.2).
- Times: `DD/MM/YYYY HH:MM:SS` (absolute); intervals `HHH:MM` (e.g. `720:00`). In a backwards
  run a source "starts" at the later time and "stops" at the earlier one, and likewise domain
  start/end times; a time grid's T0 is the chronologically **latest** time, dT stays positive
  (A01 §1.2.7.3).
- Command line `-LogFolder=<dir>` moves the Log and Error files; `Output Options: Folder` sets
  where outputs go (B01 §2, §7).
- Met: use the NAME team's declarations (`DemoScripts_v8_7_1/CommonUtilities/MetDeclarations/Use_UMG_Mk<N>_L59PTpp.txt`),
  filling the `%MetRestoreScript%`, `%MetDir%`, `%TopogDir%` placeholders, and the met
  definition file from `<NAME>/Resources/Defns/MetDefnUMG_Mk<N>_L59PTpp.txt` under `Input Files:`.
  Their met module instances are declared over the met's **whole** (global) domain.

## Domains
```
Domains:
Name, H Unbounded?, H-Coord, X Min, X Max, Y Min, Y Max, Z Unbounded?, Z-Coord, Z Max, T Unbounded?, Start Time, End Time, Max Travel Time
Dispersion Domain, No, Lat-Long, -96.329, -19.449, -64.981, 27.557, No, m agl, 20000.0, Yes, , , 720:00
```
- In Domains, X Min / X Max are the domain's **edges** (unlike grids, whose X Min is a cell
  centre); X Max < X Min is allowed for a domain crossing the dateline (B01 Domains).
- **Unbounded horizontally:** `H Unbounded? = Yes` with H-Coord **and** all four bounds left
  blank (likewise Z-Coord/Z Max blank if Z-unbounded, times blank if T-unbounded; B01): `Dispersion Domain, Yes, , , , , , No, m agl, 20000.0, Yes, , , 240:00`. Giving an
  H-Coord fails at once with "H-Coord, Set of Locations, and Location of Centre should not be
  given for unbounded domains". Unbounded runs cross the dateline and reach the poles without
  trouble.
- `Sets of Dispersion Options` names the domain as `Computational Domain` and sets the sync time
  (e.g. `00:06`). In our runs fate records appeared only at certain phases of the sync cycle
  (observation; the docs say only that they are written when a particle is removed).
- `Max # Particles` (same block) caps active particles: when reached, releases stop with
  warnings and the run carries on — leave ~20% headroom (A01 §1.5.6, §1.7.2).

## Sources
- Cuboid sources in Lat-Long × Pa: `Shape Cuboid, H-Coord Lat-Long, Z-Coord Pa`, `dX/dY` full
  widths in degrees, `dZ` full thickness in Pa, `dH-Metres? No, dZ-Metres? No`.
- `# Particles` as a rate, e.g. `30/sec`; the docs call it a lower limit, released uniformly
  through each sync interval (A01 §1.5.4.2). The number actually released by a 60-s source is not
  always 60 × rate + 1 (1,801 at 30/s, 181 at 3/s, but 600 at 10/s): count particles from the
  output, don't compute them.
- A rate on an instantaneous source (start = stop) releases nothing, and a zero maximum travel
  time makes NAME ignore the source (B03). Source `Max Age` defaults to infinity.
- At most 3,000 sources per run for footprint field requests (one request per source per output time).

## Grids and field requests
- `Horizontal Grids: Name, H-Coord, nX, nY, dX, dY, X Min, Y Min` — **X Min / Y Min are the
  centre of the first cell**, not its edge.
- **Every field request names its own source and its own horizontal grid** (`Source`, `H-Grid`
  columns). Requests in the same **output group** must use consistent grids (A01 §1.4.4.1), so
  give each grid its own output group. Different sources can be written on different grids in the same run — e.g. a small
  fine grid per release group plus shared coarse grids. Memory follows grid cells × requests
  (see performance.md), so small per-group grids are cheap and one large shared fine grid is not.
- A time-integrated footprint: `Quantity Mixing Ratio`, `T Av Or Int = Int`, one temporal grid
  with nt = 1 at the chronologically earliest end of the window, `Av Time` = window length,
  `# Av Times` = window / sync time (keep averaging steps aligned with the sync time, A01 §1.4.1.2).
  Mixing ratio needs both an H-grid and a Z-grid and can't be vertically averaged (B01).
- A blank H-Grid integrates over the whole computational domain; a blank Z-Grid integrates
  vertically (or averages over the boundary layer with `BL Average? = Yes`) (A01 §1.4.1.1).
- **A `T-Array` time grid cannot be used for field output**: the run aborts with "UNEXPECTED
  FATAL ERROR in MultiplyShortTime" (the docs say irregular time grids are not fully supported,
  A01 §1.2.6.4). Use one regular time grid (nT = 1) per time.
- **Output line length:** NAME writes each output group's fields side by side, one column per
  request. 2,720 columns in one group made the run stop **at the end** (after 8 h 21 min) with
  "FATAL ERROR: too long an output line required for output group". 2,000 columns per group works:
  split requests into output groups of at most 2,000.
- Point output on a `Locations:` set works (unstructured horizontal grid), as does a vertical
  `Z-Array`. Release-point met in one run: at most 200 time grids, 12 arrays, 86 values per array
  and 5 location sets per run (we hit these; one group per ~200 soundings).

## Particle-fate output
```
Output Requirements - Sets of Particle/Puff Information:
Output Name, Particles?, Puffs?, Fate Info?, Met?, Mass?, Plume Rise?, Dispersion Scheme?, Puff Family?, H-Coord, Z-Coord, Sync?, Output Route
Particles, Yes, No, Yes, No, No, No, No, No, Lat-Long, m agl, Yes, D
```
Writes `Particles_C1.txt`: one record per particle when it stops. `Met? = Yes` adds met columns
(incl. Pressure) to each record. Z-Coord can be `m agl` or `Pa`. See outputs.md for its quirks.

## OpenMP
`OpenMP Options: Use OpenMP?, Threads, Output Group Threads` — `Yes, 4, 1` (further columns:
Particle Threads, Particle Update Threads; A01 §1.2.7.7). Threads must equal the job's CPUs and
`OMP_NUM_THREADS`; the docs also recommend `OMP_STACKSIZE=32m`.

## Messages
A `Message Controls` block can make chosen warnings fatal (e.g. "No flow for particle/puff" =
lost particles) or stop after N messages (B01; A01 §1.2.7.8). Untested whether that changes the
exit code; we rely on checking `<input>Error.txt`.

## Met-only runs
A run that only samples met works with a single dummy source of one particle, 1 h maximum travel,
and e.g. "Pressure (Pa)" requested at 0 m agl on a `Locations:` set (one time grid per point
time): 46 points took 15 s. A pressure point below the model ground is reported at the ground (its
pressure equals the surface pressure exactly).

## Executables
`nameiii_64bit_par.exe` (text output); NetCDF field output needs
`nameiii_eccodes_NetCDF_64bit_par.exe`. Mk12 met definitions exist only in 8.7.1.
