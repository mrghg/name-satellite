"""GLIDE spec §10: checks on the one-day Mace Head runs (8 Jan 2019, 24 releases, 240 h).

Expects four run folders (``glide_runs.prepare``) under CHECK_FOLDER:
``A-fixed``, ``B-fixed-repeat`` (same input, Fixed seed), ``C-random`` (Random seed) and
``D-fixed-unpatched`` (stock NAME 8.7.1, Fixed). Writes the delivered files of run A into
``CHECK_FOLDER/delivered`` and prints each check with PASS / FAIL / INFO.

    uv run python dev/scripts/glide_tests/check10.py CHECK_FOLDER

Check 6 (identical sources in a Fixed run) was test 3 of 16 Nov 2026 (job 2296902): they differ.
"""

import sys
from pathlib import Path

import numpy as np
import xarray as xr

from name_satellite.glide_runs import spec
from name_satellite.glide_runs.deliver import deliver, footprints, read_plan
from name_satellite.release_met import read_release_met_file

MAX_AGE_TOLERANCE_HOURS = spec.SYNC_MINUTES / 60.0  # §10 check 2: 240 h plus one 300 s step


def report(label, ok, detail):
    status = "INFO" if ok is None else ("PASS" if ok else "FAIL")
    print(f"[{status}] {label}: {detail}")


def main(check_folder):
    check_folder = Path(check_folder)
    run = {name: check_folder / name for name in
           ("A-fixed", "B-fixed-repeat", "C-random", "D-fixed-unpatched")}  # fmt: skip
    out = check_folder / "delivered"
    out.mkdir(exist_ok=True)
    footprint_path, end_path = deliver([run["A-fixed"]], "MHD", "EUR-Jan", "main", out)
    ds = xr.open_dataset(footprint_path)
    end = xr.open_dataset(end_path).to_dataframe()

    # 1. Grid
    g = spec.EUR.grid
    ok = (
        ds.sizes["longitude"] == g.nx
        and ds.sizes["latitude"] == g.ny
        and np.isclose(ds.longitude[0], -29.75)
        and np.isclose(ds.longitude[-1], 39.75)
        and np.isclose(ds.latitude[0], 32.25)
        and np.isclose(ds.latitude[-1], 71.75)
    )
    report("1 grid", ok, f"{ds.sizes['longitude']} x {ds.sizes['latitude']} cells, centres "
           f"{float(ds.longitude[0])}..{float(ds.longitude[-1])}, "
           f"{float(ds.latitude[0])}..{float(ds.latitude[-1])}")  # fmt: skip

    # 2. Ages
    max_age = end["age"].max()
    report("2 oldest particle", max_age <= spec.MAX_AGE_HOURS + MAX_AGE_TOLERANCE_HOURS,
           f"{max_age:.4f} h")  # fmt: skip
    # Within a release, particle index k left the source (k - 1) / rate seconds after the end of
    # the hour (backwards run), so particles released at about xx:30 have index ~ 1,800 s x 6/s.
    middle = end[(end["particle_index"] - 10801).abs() <= 600]
    followed = middle[middle["age"] >= 239.0]
    report("2 mid-hour particles followed to 239 h", len(followed) > 0,
           f"{len(followed)} of {len(middle)} released about xx:28-xx:32 (all releases); the "
           f"rest left the domain first; oldest {middle['age'].max():.4f} h")  # fmt: skip

    # 3. Footprint totals (m2 s / mol)
    totals = ds["srr"].sum(["latitude", "longitude"]).values
    inside = ((totals >= 5) & (totals <= 17)).mean()
    report("3 footprint totals", None, f"min {totals.min():.2f}, median {np.median(totals):.2f}, "
           f"max {totals.max():.2f} m2 s/mol; {inside:.0%} within 5-17")  # fmt: skip

    # 4. Ground height at the receptor
    columns, values = read_release_met_file(run["A-fixed"] / "output" / "ReleaseMet_R1_C1.txt")
    row = dict(zip(columns["request"], values.loc["R1"], strict=True))
    report("4 receptor", None, f"model ground {row['R1_TOPO']:.2f} m asl; land fraction "
           f"{'not available' if row['R1_LAND'] < -1e30 else row['R1_LAND']}")  # fmt: skip

    # 5. Seeds
    a = footprints(run["A-fixed"], read_plan(run["A-fixed"])).values
    b = footprints(run["B-fixed-repeat"], read_plan(run["B-fixed-repeat"])).values
    c = footprints(run["C-random"], read_plan(run["C-random"])).values
    report("5 two Fixed runs identical", np.array_equal(a, b), "bit-for-bit comparison")
    report(
        "5 Fixed and Random differ",
        not np.array_equal(a, c),
        f"total differs by {np.abs(a.sum((1, 2)) - c.sum((1, 2))).max():.3g} m2 s/mol at most",
    )

    # 7. Particle counts
    n = ds["n_particles"].values
    report("7 n_particles", bool(np.all(np.abs(n - 21600) <= 5)), f"{sorted(set(n.tolist()))}")
    rows = end.groupby("time").size().values
    report("7 one end record per particle", bool(np.array_equal(rows, n)), f"{len(end)} rows")

    # 8. Patched against unpatched build
    d = footprints(run["D-fixed-unpatched"], read_plan(run["D-fixed-unpatched"])).values
    # Spec e934448: identical fate records; footprints differ by at most one unit in the 7th
    # (last printed) significant figure, i.e. a relative difference of at most 1e-6 (the 1% margin
    # allows for float rounding in the unit conversion).
    worst = np.max(np.abs(a - d) / np.maximum(np.abs(a), np.abs(d)), where=(a != 0) | (d != 0),
                   initial=0.0)  # fmt: skip
    report("8 patched = unpatched footprints", worst <= 1.01e-6,
           f"bit-identical: {np.array_equal(a, d)}; largest relative difference {worst:.2g}")  # fmt: skip

    left = ds["fraction_left_sides"].values
    top = ds["fraction_left_top"].values
    report("fractions left", None, f"sides {left.min():.2f}-{left.max():.2f}, top "
           f"{top.min():.3f}-{top.max():.3f}")  # fmt: skip


if __name__ == "__main__":
    main(sys.argv[1])
