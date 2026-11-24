# Future method: footprint file format (draft)

One netCDF-4 file per month, written by `name_satellite.future.processing`
(`<SITE>_<PRODUCT>_<YYYYMM>_future.nc`). Branch `future-method`; see `scoping.md` for the method.
Everything here is a draft for the group to test in the inversion code.

## Grids

Three output levels, each a cut-out of a **global** grid with origin (−180°, −90°):

| Level | Cell | Extent | Stored as |
|---|---|---|---|
| fine | 0.1° | ±5° around the sounding's release, snapped outward to the 2.5° lines | 125 × 125 block + offset + shape |
| medium | 0.5° | ±25° around the release, snapped outward to the 2.5° lines | 105 × 105 block + offset + shape |
| coarse | 2.5° | the globe (72 × 144) | the grid itself |

A block's cell `[j, i]` is the global cell `[offset[0] + j, offset[1] + i]` of its level, for
`j < shape[0]`, `i < shape[1]`; the rest of the block is zero padding, **not** zero sensitivity.
Because the boxes are snapped to 2.5° lines, every coarse cell is either wholly inside a finer box
or wholly outside it, and finer cells sum into coarser ones (to ~0.1%: NAME converts to mixing
ratio cell by cell).

## Variables (dimension `sounding_id` first)

| Variable | Dimensions | Units | Meaning |
|---|---|---|---|
| `srr_fine` | sounding, fine_y, fine_x | (mol/mol)/(mol m⁻² s⁻¹) | column footprint on the 0.1° block |
| `fine_offset`, `fine_shape` | sounding, yx | 1 | global 0.1° index (lat, lon) of the block's first cell; cells in the box |
| `srr_medium` | sounding, medium_y, medium_x | (mol/mol)/(mol m⁻² s⁻¹) | column footprint on the 0.5° block |
| `medium_offset`, `medium_shape` | sounding, yx | 1 | as for fine, on the global 0.5° grid |
| `srr_coarse` | sounding, coarse_lat, coarse_lon | (mol/mol)/(mol m⁻² s⁻¹) | column footprint on the global 2.5° grid |
| `end_positions` | sounding, end_height, coarse_lat, coarse_lon | 1 | where the particles were at the end of the back-run: column-weighted fraction per 2.5° × 1 km bin (0–20 km agl); the last height = left through the 20 km top. Sums to the sounding's column weight |
| `end_time` | sounding | time | when particles reached the maximum travel time (240 h before the release) |
| `top_exit_mean_age` | sounding | hours | mean travel time of the particles that left through the top |
| `column_weight` | sounding, release_level | 1 | pressure weight × averaging kernel of each modelled level |
| `time`, `lat`, `lon` | sounding | | the sounding (from the observation file) |
| `release_group`, `release_lat`, `release_lon`, `release_time` | sounding | | the shared release (D20) |
| `fine_lat_global`, `fine_lon_global`, `medium_lat_global`, `medium_lon_global` | | degrees | cell centres of the global 0.1° and 0.5° grids |

Footprints are the column footprint of each sounding, Σ over modelled levels of column weight ×
the level's footprint (one per sounding, D26). Using them:

- **Enhancement**: Σ over the cells of the inversion's grid of footprint × flux, taking for each flux
  cell the finest level the sounding has there; where that level is coarser than the flux cell,
  spread its value evenly (`scoping.md`, "Using the nested grids in an inversion").
- **Boundary term**: Σ over `end_positions` bins of the fraction × the CH₄ mole fraction of a global
  3D field in that bin at `end_time` (top exits: at the 20 km top, `top_exit_mean_age` hours before
  the release).
- `name_satellite.future.processing.sounding_footprint(ds, sounding_id, level)` returns one
  sounding's footprint on one level with its latitude and longitude, as an example of placing the
  blocks.

Storage: float32 rounded to 10 significant bits (BitRound, unbiased, 0.05% relative precision),
zlib 4 + shuffle, one sounding per chunk.
