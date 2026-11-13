"""How much would the soundings of a release group lose by sharing one footprint? (D20, D23)

Each sounding's column footprint is Σ_j w_sj F_gj: its own weights w_sj (pressure weight ×
averaging kernel) on its group's per-level footprints F_gj (processing.process_day). The file can
point several soundings at one footprint (``footprint_index``, D23). The natural shared footprint
uses the group-mean weights, Σ_j w̄_j F_gj, which by linearity is the mean of the members' own
columns. This script compares each member's own column with that shared one, from the unrounded
``process_day`` output of one day:

* "mean": the shared footprint as it is;
* "mean, rescaled": the shared footprint scaled by Σ_j w_sj / Σ_j w̄_j, so that its total
  particle fraction (the boundary weight) equals the sounding's own column weight.

Reported per member of a multi-sounding group:

* srr total and cell-by-cell (L1) differences, relative to the own footprint, and as XCH4 for an
  illustrative uniform flux (FLUX_MOL_M2_S, roughly South America's CH4 emissions spread evenly);
* the boundary term: Σ over faces, top and inside of the difference, times a uniform boundary
  mole fraction BOUNDARY_PPB (the error in the modelled background if the inversion used the
  shared sensitivities with no correction).

Usage (sci node; one day reads ~0.5 GB of NAME fields):
    uv run python dev/scripts/measure_grouping_error.py DAY_FOLDER
"""

import sys
from pathlib import Path

import numpy as np

from name_satellite.processing import FACES, process_day

FLUX_MOL_M2_S = 1e-8  # ~100 Tg CH4/yr over ~1.8e7 km2: illustrative only
BOUNDARY_PPB = 1850.0
PPB = 1e9


def boundary_total(ds, i):
    """Total particle fraction (all faces, top and inside) of footprint i."""
    total = float(ds["particle_fraction_inside"][i]) + float(ds["particle_locations_top"][i].sum())
    for face in FACES:
        total += float(ds[f"particle_locations_{face}"][i].sum())
    return total


def boundary_l1(ds, i, shared, scale):
    """Σ |scale × shared − own| over every boundary array of footprint i."""
    names = [f"particle_locations_{f}" for f in FACES] + ["particle_locations_top"]
    names.append("particle_fraction_inside")
    return sum(float(np.abs(scale * shared[name] - ds[name][i].values).sum()) for name in names)


def main(day_folder):
    ds = process_day(Path(day_folder), allow_test=False)
    groups = ds["release_group"].values
    weight_sum = ds["column_weight"].sum("release_level").values
    rows = []
    for group in np.unique(groups):
        members = np.flatnonzero(groups == group)
        if len(members) < 2:
            continue
        names = ["srr", "particle_fraction_inside", "particle_locations_top"]
        names += [f"particle_locations_{f}" for f in FACES]
        shared = {n: ds[n].values[members].mean(axis=0) for n in names}
        shared_weight = weight_sum[members].mean()
        for i in members:
            own = ds["srr"].values[i]
            for variant, scale in (("mean", 1.0), ("rescaled", weight_sum[i] / shared_weight)):
                fp = scale * shared["srr"]
                rows.append(
                    {
                        "group": group,
                        "n": len(members),
                        "variant": variant,
                        "weight_diff": shared_weight * scale - weight_sum[i],
                        "total_rel": fp.sum() / own.sum() - 1.0,
                        "l1_rel": np.abs(fp - own).sum() / own.sum(),
                        "enh_ppb": (fp.sum() - own.sum()) * FLUX_MOL_M2_S * PPB,
                        "own_enh_ppb": own.sum() * FLUX_MOL_M2_S * PPB,
                        "bnd_ppb": (scale * sum(float(shared[n].sum()) for n in names[1:])
                                    - boundary_total(ds, i)) * BOUNDARY_PPB,
                        "bnd_l1": boundary_l1(ds, i, shared, scale) / boundary_total(ds, i),
                    }
                )  # fmt: skip
    print(f"{len(ds.sounding_id)} soundings, {len(np.unique(groups))} groups; "
          f"{sum(r['variant'] == 'mean' for r in rows)} soundings in multi-member groups")  # fmt: skip
    members = [r for r in rows if r["variant"] == "mean"]
    weight_diff = np.abs([r["weight_diff"] for r in members])
    print(
        f"column weight, shared - own: median |Δ| {np.median(weight_diff):.2e}, "
        f"max {weight_diff.max():.2e}"
    )
    for variant in ("mean", "rescaled"):
        sel = [r for r in rows if r["variant"] == variant]
        print(f"\n{variant}:")
        for key, label, fmt in (
            ("total_rel", "srr total, relative", "{:.2e}"),
            ("l1_rel", "srr cell-by-cell L1, relative", "{:.2e}"),
            ("enh_ppb", "enhancement error (ppb, uniform flux)", "{:.4f}"),
            ("bnd_l1", "boundary arrays L1, relative", "{:.2e}"),
            ("bnd_ppb", "boundary term error (ppb, uniform 1850)", "{:.3f}"),
        ):
            v = np.abs([r[key] for r in sel])
            print(f"  {label:42s} median |.| " + fmt.format(np.median(v))
                  + "   90% " + fmt.format(np.percentile(v, 90)) + "   max " + fmt.format(v.max()))  # fmt: skip
    own = np.array([r["own_enh_ppb"] for r in members])
    print(
        f"\nown enhancement at that flux: median {np.median(own):.1f} ppb "
        f"({own.min():.1f}-{own.max():.1f})"
    )


if __name__ == "__main__":
    main(sys.argv[1])
