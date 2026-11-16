# Hand-off: NAME reference runs for GLIDE

The spec is `~/code/glide/docs/reference-runs.md` (read it in full; its numbers are fixed). It asks
for NAME back-run footprints from 7 point receptors (one a column) on two 0.25° grids (EUR, NAM),
hourly releases over two 14-day blocks in 2019, plus repeat-seed and single-release runs.

We compared three options: enw, a bespoke tool, and this repo. **Build it here, on a branch off
`main`** (not `future-method`: the spec needs a bounded box with particles removed at its edge).
Most of the pipeline already does what the spec needs:

- one source and one integrated 0–40 m mixing-ratio field per release;
- a box domain with a fixed per-particle `Max Travel Time` (set it to 240 h), with a 20 km top (D21);
- Mk10 met (covers Dec 2018 – Jul 2019), met restore, LOTUS jobs;
- particle-fate records from the patched NAME, so `fraction_left_sides`, `fraction_left_top` and
  the §8 end records come from `classify_fates`/`release_fates`;
- model surface pressure at a point (`surface_pressure.py`) and pressure layers (`releases.py`)
  for the Lamont column;
- the conversion to (mol/mol)/(mol m⁻² s⁻¹) (`footprints.py`).

## Changes needed

1. **Driver.** A new entry point to use in place of the GOSAT-specific `prepare.py`. It builds the
   release plan from the receptor table (§3) × hourly release times (§4), then reuses
   `write_name_input`, `lotus.py` and `met.py`. Keep it separate from the satellite path.
2. **Point sources in m agl.** `sources_and_requests_block` hard-codes cuboids in `Pa`. Allow the
   Z-coord to be `m agl` with dX = dY = dZ = 0. Use 20,000 particles per release, spread through
   the hour (100,000 for Lamont).
3. **Lamont column.** From the model surface pressure up to 100 hPa (§5). First check whether a
   cuboid in `Pa` spreads particles evenly in pressure or evenly in height. If it is even in
   pressure, one source is enough. If not, use 30 stacked layers. Record the method and the layer
   boundaries in metres.
4. **Convection.** `dispersion_options_block` has `Deep Convection? No`. The spec (§7) wants it on,
   as in the group's standard setup. Make it a `RunSettings` field.
5. **Grids and domains.** Add EUR and NAM exactly as in spec §2 (279 × 159 and 259 × 143 cells,
   first cell centres given there).
6. **Single-release runs (§9).** Particle positions at 6, 24, 72 and 120 h of travel. Add a
   particle-output request on its own temporal grid, next to the existing fate output.
7. **Output writer.** Write one file per receptor, block and run, in the §8 layout (`srr`,
   `release_lon/lat/height`, Lamont pressure bounds, `n_particles`, `fraction_left_*`), plus an
   `_end` file, using the §8 file names. Reuse the readers and the fate code. Don't reuse the
   column-combining part of `processing.py`.
8. **Run record (§8).** Record ground height and surface type at each receptor (check whether
   release-met can give land/sea type), the Lamont surface pressures, and the domain top and
   what happens to particles there.

## Before running

- **Random seeds:** NAME has no numeric seed. The GLIDE side has been told, and the spec owner
  will decide how §9's repeat-seed runs are done. Don't improvise a scheme.
- **Follow spec §10:** run one day at MHD first, and send the run record before starting the blocks.
- **Licence:** NAME output and Met Office fields never go into the GLIDE repository.
