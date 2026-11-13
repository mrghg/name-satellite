"""Figures for the group slides on the first full test day of the new workflow (November 2026).

Reads the processed footprint file of the new 1 Jan 2013 test day (processing.py, D23) and the
legacy production file for January 2013, and writes two PNGs:

    f7_newrun_day.png       sum of the day's 46 column footprints, with the sounding positions
    f8_exit_shares.png      share of each sounding's side exits through each face: new run vs the
                            nearest legacy production entry (1:1 line), one panel per face

The style follows make_slide_figures.py. matplotlib and cartopy are not project dependencies:

    uv run --with matplotlib --with cartopy python dev/scripts/figures/make_newrun_figures.py \
        NEW_FILE OUTPUT_FOLDER
"""

import sys
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from name_satellite.validation.refdata import ref_path

FP_UNITS = "(mol/mol) / (mol m⁻² s⁻¹)"
FACES = {"n": "North", "e": "East", "s": "South", "w": "West"}
INK = "#1f2430"
MUTED = "#52514e"
SURFACE = "#fbfaf7"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
SEQUENTIAL = mcolors.LinearSegmentedColormap.from_list(
    "sequential", ["#e6f0fc", "#9ec5f4", "#5598e7", "#256abf", "#104281", "#0d366b"]
)
SEQUENTIAL.set_under(SURFACE)
plt.rcParams.update(
    {
        "font.size": 18,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
    }
)


def side_exit_shares(fluxes):
    """{face: share of the side exits} from {face: per-sounding total}."""
    total = sum(fluxes.values())
    return {face: value / total for face, value in fluxes.items()}


def figure_day(out, new):
    """F7: sum of the day's column footprints, with the sounding positions."""
    srr = new["srr"].isel(footprint=new["footprint_index"])
    lat, lon = new["latitude"].values, new["longitude"].values
    total = srr.sum("sounding_id").values
    vmax = float(total.max())
    fig = plt.figure(figsize=(9, 11))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.8, edgecolor=MUTED)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.5, edgecolor="#9a9890")
    mesh = ax.pcolormesh(
        lon, lat, total, norm=mcolors.LogNorm(vmax * 1e-4, vmax), cmap=SEQUENTIAL,
        transform=ccrs.PlateCarree(),
    )  # fmt: skip
    ax.scatter(
        new["release_lon"], new["release_lat"], s=36, facecolor="none", edgecolor=ORANGE,
        linewidth=1.5, transform=ccrs.PlateCarree(), label="Soundings",
    )  # fmt: skip
    ax.legend(loc="lower left", frameon=False)
    colorbar = fig.colorbar(mesh, ax=ax, shrink=0.6, pad=0.02, extend="min")
    colorbar.set_label(f"Sum of the day’s column footprints,\n{FP_UNITS}")
    fig.savefig(out / "f7_newrun_day.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def figure_exit_shares(out, new, legacy):
    """F8: share of side exits by face, new vs nearest legacy entry (same day)."""
    order = new["footprint_index"].values
    new_flux = {
        f: new[f"particle_locations_{f}"].sum(("height", "longitude" if f in "ns" else "latitude"))
        .values[order]
        for f in FACES
    }  # fmt: skip
    old_flux = {
        f: legacy[f"particle_locations_{f}"].sum(("height", "lon" if f in "ns" else "lat")).values
        for f in FACES
    }
    new_share, old_share = side_exit_shares(new_flux), side_exit_shares(old_flux)

    # Nearest legacy entry: legacy soundings were binned, so several new ones map to one.
    seconds = (new["time"].values[:, None] - legacy["time"].values[None, :]) / np.timedelta64(
        1, "s"
    )
    distance = np.hypot(
        new["release_lat"].values[:, None] - legacy["release_lat"].values[None, :],
        new["release_lon"].values[:, None] - legacy["release_lon"].values[None, :],
    )
    nearest = np.argmin(np.abs(seconds) / 60 + 10 * distance, axis=1)
    if (np.abs(seconds[np.arange(len(nearest)), nearest]) > 600).any():
        raise ValueError("A new sounding has no legacy entry within 10 minutes")

    fig, axes = plt.subplots(1, 4, figsize=(20, 5.6), sharex=True, sharey=True)
    for ax, (face, name) in zip(axes, FACES.items(), strict=True):
        x, y = old_share[face][nearest], new_share[face]
        ax.plot([0, 1], [0, 1], color="#9a9890", linewidth=1.5, zorder=1)
        ax.scatter(x, y, s=70, color=BLUE, edgecolor=SURFACE, linewidth=1.2, zorder=2)
        ax.set_title(name, loc="left", fontsize=20)
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.03, 1.03)
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks([0, 0.5, 1])
        ax.set_aspect("equal")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        diff = np.abs(y - x)
        ax.text(0.97, 0.05, f"median |Δ| {np.median(diff):.3f}", ha="right", fontsize=15,
                color=MUTED, transform=ax.transAxes)  # fmt: skip
    axes[0].set_ylabel("New run")
    fig.supxlabel("Legacy production (nearest entry)", fontsize=18, y=0.0)
    fig.savefig(out / "f8_exit_shares.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def main(new_file, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    new = xr.open_dataset(new_file).load()
    legacy = xr.open_dataset(ref_path("production_201301")).sel(time="2013-01-01").load()
    figure_day(out, new)
    figure_exit_shares(out, new, legacy)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
