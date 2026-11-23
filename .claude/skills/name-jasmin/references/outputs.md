# NAME 8.7.1 output files

Readers: `~/code/name-satellite/src/name_satellite/name_output.py` (`read_fields`,
`read_name_fate_file`, `drop_repeated_fate_records`).

## Files next to the input
- `<input>Log.txt`: ends with "Run completed successfully at ..." on success.
- `<input>Error.txt`: must be empty. **NAME's exit code is 0 even when errors were written here.**
- `<input>Timing.txt`: timer summary (with `Timer Options: Use Timers? Yes`); the first block is
  the total run time in seconds.
- `sacct` TotalCPU shows 0 under mpirun; MaxRSS is sampled every ~30 s, so short runs show only
  mpirun's few MB (see jasmin-lotus.md for measuring NAME's own peak memory).

## Fields files (`<OutputGroup>_C1_T1_<YYYYMMDDHHMM>.txt`, Output Format `IA`)
- Header of `Key: value` lines (run name, grid origin/size/resolution, number of field columns),
  then `Fields:` and one header row per piece of column metadata (species category, request name,
  quantity, species, units e.g. "ppm s", **source name**, averaging...), then `X Index, Y Index,
  X, Y, value...` rows.
- **Only grid cells where some column is non-zero are written**; restore the rest as zeros.
  Indices are 1-based. "X/Y grid origin" = the input's X Min / Y Min = centre of cell 1.
- Coordinates printed to 5 decimals; lines end with a trailing comma.
- Match columns to releases by the source-name header row, never by column position.
- Mixing ratio integrated over the window (ppm s) → footprint in (mol/mol)/(mol m⁻² s⁻¹) with a
  1 g/s source of molecular weight 1: field × 1e-6 × cell area (m²) / release duration (s).
- **How air density enters `Mixing Ratio`** (test of 23 Nov 2026, job 6142203; Mk10, Mace Head
  case-1 release, 24 h, 0–40 m layer; instantaneous `Air Concentration`, `Density`,
  `Temperature (K)`, `Pressure (Pa)` and `Mixing Ratio` on 8 × 8 cells at all 301 sync times,
  plus time-integrated requests):
  - NAME divides the cell's air concentration by the air density **at the centre of the output
    grid box** (documented: A01 output quantities; C03 §2.2), **at each sampling time, before
    time integration**: the integrated mixing ratio equals Σ (C/ρ) Δt × 28.966 / M from the
    instantaneous values to 4e-7, while (∫C dt) / mean ρ differs by 0.15–0.4%.
  - The density is p / (287.05 T) with the actual temperature, to 1e-6: **no humidity term**,
    so the "dry air" mixing ratio is in effect per mole of moist air. Molar mass of air 28.966.
  - So: one density per output cell and sampling time, never per particle. Against a model that
    weights each particle by its own 1/ρ, expect differences of a few percent in steep terrain
    (density at the cell-centre ground height) and from the dry/moist-air definition (≈ e/p, up
    to 1–2.5% near the surface in humid summer air).
  - A `Density` request must have `Sync? = Yes`, or the run stops at once with a fatal error
    saying so (also when time-integrated).
- Output Format letters (B01; A01 Table 1.12): I = include grid indices, A = align columns,
  **Z = also write all-zero grid points** (default: omit them), 2 = legacy NAME II layout
  (deprecated). `Masking Threshold` can trim small values further.
- Sizes measured (48 sources, 0–40 m layer, time-integrated): 0.1° ±5° box per group 0.5–3 MB;
  0.5° over South America 13 MB (240 h) / 44 MB (720 h); 2.5° globe 3–15 MB.

## Particle-fate file (`Particles_C1.txt`, "Fate Info? = Yes")
- After a header, `Source,Travel Time (hr),X (Lat-Long),Y (Lat-Long),Z (<units>)` and one record
  per particle when it stops: leaves the domain (side or top) or reaches the maximum travel time.
  With `Met? = Yes`, met columns follow and lines end with a comma.
- **Stock 8.7.1 fate records carry no particle identifier and no stop reason** (as B01 lists
  them). Stock NAME's ordinary trajectory output (not Fate Info) does carry a per-source particle
  index (A01 §1.4.2). Particles are also removed at the domain end time or the source's Max Age
  (A01 §1.2.6.5). Fate output can be limited with Source and First/Last Particle columns (B01).
- **Stock 8.7.1 writes some particles twice** (~2–3% of records in 30-day runs, almost all at
  domain exits; a known issue reported to the Met Office). Counts per source therefore can't be
  used to count particles released.
- The ACRG-patched build adds `Particle Index` (1..n per source) and `Fate Code` after Z. Fate
  codes are bit flags: 8 = left the domain, 16 = reached the maximum travel time; 24 = both (left
  on its last step).
- With the patched build, repeated records are always same source, same index, both "left the
  domain", one particle time step apart (10 s in a small test; **180 s** in a 352-source 30-day
  run). The first record is the exit (always outside the domain). Keep each particle's first
  record (allow gaps up to the sync time).
- A few "left the domain" records sit exactly on an edge as printed (beyond it by less than the
  5-decimal precision): treat as exits through that edge.
- In unbounded runs (no side boundaries) every particle has one record at the maximum travel time
  (fate 16) except those leaving through the top (fate 8): these are the 3D end positions.
- Size: ~545,000 records ≈ 157 MB for 352 sources with met columns.

## Point output (`<group>_C1.txt`, on a `Locations:` set)
Fields-file layout with 15 header rows; preliminary columns `T Index, Location Index, T,
Location`; one row per location; 7 significant figures.

## Particle positions on a time grid (`<Output Name>_C1.txt`, particle output with a `T-Grid`)
Test of 16 Nov 2026 (job 2295092; GLIDE smoke test, 3 one-hour releases at 6/s, 1 h travel,
backwards, hourly grid of 4 times from T0 = end of the last release):
- A `Sets of Particle/Puff Information` row with `Fate Info? = No` and a `T-Grid` writes every
  live particle at each grid time. Columns: `Puff?, PP Index` (per-source particle number, the
  same index as the patched fate records), `Source Index, Source, Release Time, Time, Travel Time`
  (**seconds**), `X, Y, Z` (in the requested Z-Coord), turbulent velocities and puff spreads.
  10 header lines before the column heading.
- In a backwards run the grid counts back from T0 (as for fields). Nothing is written at T0
  itself when no particle exists yet.
- A particle that reaches the maximum travel time stops moving at exactly that age (fate record
  age = 1.000 h for all 64,803 particles) but stays listed, frozen at that age, until the end of
  that sync step; it is removed after. So output at a grid time can include particles held at the
  maximum age (1,800 = one 5-min sync step at 6/s).
- A one-hour source at `6/sec` released **21,601** particles (rate × 3,600 + 1), in both the
  fate and position outputs.

## Met quantities at a point: availability with UM Mk10 met
Same test: `Topography` (m) and `Pressure (Pa)` / `Boundary Layer Depth` at 0 m agl were
written (Mace Head: 1.32 m, 103,462.5 Pa, 134 m). **`Land Fraction` came back as -3.402823E+38**
(the missing-value marker): the Mk10 met does not supply it, so NAME cannot report a point's land
or sea type from these fields.
