"""Figures for the group slides on the stage 1 and 2 validation (November 2026).

Reads the processed footprint files of the stage 1 runs (A, B, B2) and the stage 2 daily runs
(written by process_stage1.py and process_stage2.py) and the production file, and writes four PNGs:

    f1_footprint_compare.png  one sounding: production vs rebuild (seed 1) vs rebuild (seed 2)
    f2_totals_scatter.png     column totals: rebuild vs production (1-7 Jan), and seed 1 vs seed 2
    f3_x18_window.png         X18: change from correcting the integration window vs release delay
    f4_day_footprints.png     sum of the 1 Jan column footprints from run A (corrected window)

matplotlib and cartopy are not project dependencies, so run with:

    uv run --with matplotlib --with cartopy python dev/scripts/make_slide_figures.py OUTPUT_FOLDER
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
import pandas as pd
import xarray as xr

from name_satellite.refdata import ref_path
from name_satellite.validation import near_field_ratio

SCRATCH = Path("/work/scratch-pw5/mattrigby/name-satellite")
STAGE1 = SCRATCH / "stage1_20261109"
STAGE2 = SCRATCH / "stage2_20261109"
SOUNDING = 13  # 1 Jan sounding shown in F1 (large footprint, released 1.66 h after the first)

# Colours: one blue sequential ramp for magnitudes; blue and orange for the two series.
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


def read_fp(path):
    """Footprints (time, lat, lon) from a processed (legacy-layout) footprint file."""
    return xr.open_dataset(path)["fp"].transpose("time", "lat", "lon").load()


def map_axes(fig, position, extent):
    ax = fig.add_subplot(*position, projection=ccrs.PlateCarree())
    ax.set_extent(extent, ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.8, edgecolor=MUTED)
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.5, edgecolor="#9a9890")
    return ax


def tidy(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def figure_footprints(out, production, b, b2):
    """F1: one sounding in production, rebuild seed 1 and rebuild seed 2 (log scale)."""
    s = SOUNDING
    lat, lon = b["lat"].values, b["lon"].values
    lon0 = float(lon[np.argmax(b[s].max("lat").values)])
    lat0 = float(lat[np.argmax(b[s].max("lon").values)])
    extent = [lon0 - 9, lon0 + 9, lat0 - 7, lat0 + 7]
    panels = [
        (production[s], "2019 production\n(NAME III v7.2)"),
        (b[s], "Rebuild, legacy settings\n(NAME 8.7.1, seed 1)"),
        (b2[s], "Same, different seed\n(NAME 8.7.1, seed 2)"),
    ]
    vmax = float(max(p.max() for p, _ in panels))
    norm = mcolors.LogNorm(vmin=vmax * 1e-3, vmax=vmax)
    fig = plt.figure(figsize=(19, 7.2))
    for i, (data, title) in enumerate(panels):
        ax = map_axes(fig, (1, 3, i + 1), extent)
        mesh = ax.pcolormesh(
            lon, lat, data.values, norm=norm, cmap=SEQUENTIAL, transform=ccrs.PlateCarree()
        )
        ax.set_title(title, fontsize=20, pad=10)
    colorbar = fig.colorbar(mesh, ax=fig.axes, shrink=0.8, pad=0.02, extend="min")
    colorbar.set_label("Column footprint, (mol/mol) / (mol m⁻² s⁻¹)")
    fig.savefig(out / "f1_footprint_compare.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def figure_totals(out, production_all, production_0101, b, b2):
    """F2: column totals, rebuild vs production (all days available) and seed 1 vs seed 2."""
    new_totals = [b.sum(("lat", "lon")).values]
    production_totals = [production_0101.sum(("lat", "lon")).values]
    last_day = 1
    for d in range(2, 8):
        day = f"201301{d:02d}"
        path = STAGE2 / f"day_{day}" / f"fp_{day}_legacy_processing.nc"
        if path.exists():
            daily = read_fp(path)
            new_totals.append(daily.sum(("lat", "lon")).values)
            production_day = production_all.sel(time=daily["time"].values)
            production_totals.append(production_day.sum(("lat", "lon")).values)
            last_day = d
    new_totals = np.concatenate(new_totals)
    production_totals = np.concatenate(production_totals)
    b_totals = b.sum(("lat", "lon")).values
    b2_totals = b2.sum(("lat", "lon")).values
    scatter_rebuild = 100 * np.log(new_totals / production_totals).std()
    scatter_seeds = 100 * np.log(b_totals / b2_totals).std()
    print(f"F2: {len(new_totals)} soundings, rebuild/production median "
          f"{np.median(new_totals / production_totals):.4f}, scatter {scatter_rebuild:.1f}%; "
          f"seeds scatter {scatter_seeds:.1f}%")  # fmt: skip

    everything = np.concatenate([new_totals, production_totals, b2_totals])
    limits = [everything.min() * 0.7, everything.max() * 1.4]
    panels = [
        (production_totals, new_totals, "2019 production", "Rebuild (NAME 8.7.1)", BLUE,
         f"Rebuild vs production, 1–{last_day} Jan 2013\n"
         f"{len(new_totals)} soundings, scatter {scatter_rebuild:.1f}%"),
        (b2_totals, b_totals, "Rebuild, seed 2", "Rebuild, seed 1", ORANGE,
         f"Seed 1 vs seed 2 (noise), 1 Jan 2013\n25 soundings, scatter {scatter_seeds:.1f}%"),
    ]  # fmt: skip
    fig, axes = plt.subplots(1, 2, figsize=(15, 7.2))
    for ax, (x, y, x_label, y_label, colour, title) in zip(axes, panels, strict=True):
        ax.plot(limits, limits, color="#9a9890", linewidth=1.5, zorder=1)
        ax.scatter(x, y, s=60, color=colour, edgecolor=SURFACE, linewidth=1.2, zorder=2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.set_aspect("equal")
        ax.set_xlabel(f"{x_label}: column total")
        ax.set_ylabel(f"{y_label}: column total")
        ax.set_title(title, fontsize=20, pad=10)
        tidy(ax)
    fig.tight_layout()
    fig.savefig(out / "f2_totals_scatter.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def figure_window(out, a, b, positions):
    """F3 (X18): runs A (corrected window) vs B (legacy window), same seed, by release delay."""
    times = pd.to_datetime(a["time"].values)
    hours = (times - times.min()).total_seconds() / 3600
    total_change = 100 * (a.sum(("lat", "lon")).values / b.sum(("lat", "lon")).values - 1)
    near = near_field_ratio(
        a.values, b.values, a["lat"].values, a["lon"].values,
        positions["release_lat"].values, positions["release_lon"].values,
    )  # fmt: skip
    near_change = 100 * (near - 1)
    print(f"F3: total change {total_change.min():.2f} to {total_change.max():.2f}%, "
          f"near field {near_change.min():.1f} to {near_change.max():.1f}%")  # fmt: skip

    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), sharex=True)
    panels = [
        (total_change, "Whole footprint", BLUE),
        (near_change, "Within 2° of the sounding", ORANGE),
    ]
    for ax, (change, title, colour) in zip(axes, panels, strict=True):
        ax.axhline(0, color="#9a9890", linewidth=1.5, zorder=1)
        ax.scatter(hours, change, s=80, color=colour, edgecolor=SURFACE, linewidth=1.2, zorder=2)
        ax.set_title(title, fontsize=22, pad=10)
        ax.set_xlabel("Release time after the day’s first sounding (h)")
        ax.set_ylabel("Change, corrected vs legacy window (%)")
        tidy(ax)
    fig.tight_layout()
    fig.savefig(out / "f3_x18_window.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def figure_day(out, a, positions):
    """F4: sum of the day's column footprints from run A, with the sounding positions."""
    lat, lon = a["lat"].values, a["lon"].values
    total = a.sum("time").values
    vmax = float(total.max())
    fig = plt.figure(figsize=(9, 11))
    ax = map_axes(fig, (1, 1, 1), [lon.min(), lon.max(), lat.min(), lat.max()])
    mesh = ax.pcolormesh(
        lon, lat, total, norm=mcolors.LogNorm(vmax * 1e-4, vmax), cmap=SEQUENTIAL,
        transform=ccrs.PlateCarree(),
    )  # fmt: skip
    ax.scatter(
        positions["release_lon"], positions["release_lat"], s=36, facecolor="none",
        edgecolor=ORANGE, linewidth=1.5, transform=ccrs.PlateCarree(), label="Soundings",
    )  # fmt: skip
    ax.legend(loc="lower left", frameon=False)
    colorbar = fig.colorbar(mesh, ax=ax, shrink=0.6, pad=0.02, extend="min")
    colorbar.set_label("Sum of the day’s column footprints")
    fig.savefig(out / "f4_day_footprints.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    a_file = STAGE1 / "run_A" / "fp_20130101_legacy_processing.nc"
    a = read_fp(a_file)
    b = read_fp(STAGE1 / "run_B" / "fp_20130101_legacy_processing.nc")
    b2 = read_fp(STAGE1 / "run_B2" / "fp_20130101_legacy_processing.nc")
    positions = xr.open_dataset(a_file)[["release_lat", "release_lon"]].load()
    production_all = xr.open_dataset(ref_path("production_201301"))["fp"]
    production_all = production_all.transpose("time", "lat", "lon")
    production_0101 = production_all.sel(time=b["time"].values).load()

    figure_footprints(out, production_0101, b, b2)
    figure_totals(out, production_all, production_0101, b, b2)
    figure_window(out, a, b, positions)
    figure_day(out, a, positions)
    print(f"Figures written to {out}")


if __name__ == "__main__":
    main(sys.argv[1])
