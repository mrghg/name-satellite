"""Reproduce the NAME source definitions of the 2013-01-01 legacy GOSAT-BRAZIL run.

Re-implements acrg.satellite.common.midpoint_bounds (set_as_edges=False) and compares it with the
425 `SourceID_*` lines in the legacy NAME input file.

Finding (2026-11-06): levels 3-17 match using the GOSAT pressure levels alone; levels 1-2 match only
when p0 is replaced by a surface pressure inferred from the input (the legacy NAME surface-pressure
run, use_name_pressure=True).

Run:  uv run python dev/scripts/legacy_validation/reproduce_legacy_sources.py
"""

import numpy as np
import xarray as xr

from name_satellite.validation.refdata import ref_path

N_LEVELS = 17


def midpoint_layers(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Legacy release layers (centre, thickness) from level pressures p[n_soundings, n_levels]."""
    d = p[:, :-1] - p[:, 1:]
    dp = np.zeros(p.shape)
    centre = np.zeros(p.shape)
    dp[:, 0] = d[:, 0] / 2
    centre[:, 0] = p[:, 0] - dp[:, 0] / 2
    dp[:, 1:-1] = d[:, :-1] / 2 + d[:, 1:] / 2
    dp[:, -1] = d[:, -1]
    below = p[:, -1] - dp[:, -1] / 2 < 0
    dp[below, -1] = (p[:, -1] + d[:, -1] / 2)[below]
    centre[:, 1:] = (p[:, 1:] + d / 2) - dp[:, 1:] / 2
    return centre, dp


def read_sources(path):
    rows = [
        [c.strip() for c in line.split(",")]
        for line in path.read_text().splitlines()
        if line.startswith("SourceID_")
    ]
    cols = {"X": 12, "Y": 13, "Z": 14, "dX": 15, "dY": 16, "dZ": 17}
    return {k: np.array([float(r[i]) for r in rows]) for k, i in cols.items()}


def main():
    src = read_sources(ref_path("mr2019_input"))
    obs = xr.open_dataset(ref_path("obs"))
    n = obs.sizes["time"]
    assert src["Z"].size == n * N_LEVELS, "unexpected number of sources"
    z = src["Z"].reshape(n, N_LEVELS)
    dz = src["dZ"].reshape(n, N_LEVELS)

    p = obs.pressure_levels.values.astype(float) * 100.0  # hPa -> Pa

    # Surface pressure implied by the legacy bottom layer: centre = ps - dz0/2
    ps_name = z[:, 0] + dz[:, 0] / 2
    print("Inferred NAME surface pressure minus GOSAT p0 (hPa):")
    print(np.round((ps_name - p[:, 0]) / 100, 2))

    for label, p0 in [("GOSAT p0", p[:, 0]), ("NAME ps", ps_name)]:
        pp = p.copy()
        pp[:, 0] = p0
        centre, thick = midpoint_layers(pp)
        rel_z = np.abs(z - centre[:, :N_LEVELS]) / centre[:, :N_LEVELS]
        rel_dz = np.abs(dz - thick[:, :N_LEVELS]) / thick[:, :N_LEVELS]
        print(f"{label:9s}: max rel diff Z = {rel_z.max():.2e}, dZ = {rel_dz.max():.2e}")

    for name, var in [("X", "lon"), ("Y", "lat")]:
        diff = np.abs(src[name] - np.repeat(obs[var].values, N_LEVELS)).max()
        print(f"{name} max abs diff vs obs: {diff:.2e}")


if __name__ == "__main__":
    main()
