"""Results and figures of the future-method test runs (dev/future/prepare_unbounded_test.py).

Part of the scoping work on the ``future-method`` branch (docs/future/scoping.md).

Reads the smoke runs (1 h travel) and the full runs (720 h, 240 h) of the three test release
groups, prints a summary and writes PNGs for the group slides:

    fu1_smoke_cost.png        peak memory and run time of each smoke variant (1 h travel)
    fu2_full_cost.png         run time and peak memory of the full runs
    fu3_end_positions.png     where the particles are when they stop: today's runs (domain exits)
                              vs the unbounded runs at 240 h and 720 h
    fu4_nested_footprint.png  one release group's column footprint on the nested grids (0.1° per
                              group, 0.5° region, 2.5° globe) as sensitivity per unit area
    fu5_region_vs_outside.png share of each group's column footprint inside the region and
                              outside it, and today's bounded footprint against the unbounded one
    fu6_by_distance.png       each group's column footprint by distance from its release (inside
                              the region), today vs 240 h vs 720 h
    summary.json              the numbers behind the figures

Column footprints use the release group's own sounding (its ID is the group ID): footprint of
each level × pressure weight × averaging kernel, summed over the 16 modelled levels, as
processing.process_day. Footprints from different grids are compared as sensitivity per unit area
(footprint / cell area), which does not depend on the cell size.

matplotlib and cartopy are not project dependencies:

    uv run --with matplotlib --with cartopy python dev/future/analyse_unbounded_test.py \\
        SMOKE_FOLDER FULL_FOLDER OBS_FILE OUTPUT_FOLDER
"""

import json
import re
import subprocess
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

from name_satellite.footprints import footprint_from_mixing_ratio, grid_cell_area_m2
from name_satellite.name_output import drop_repeated_fate_records, read_fields, read_name_fate_file

N_LEVELS = 16
RELEASE_SECONDS = 60.0
SYNC_SECONDS = 360.0
# Today's dispersion domain (domains.SOUTHAMERICA) and the region of the test's medium grid.
DOMAIN = {"x_min": -96.329, "x_max": -19.449, "y_min": -64.981, "y_max": 27.557}
REGION = {"x_min": -97.0, "x_max": -19.0, "y_min": -65.0, "y_max": 28.0}
GROUP_NAMES = {1: "NE coast of Brazil", 2: "Amazon", 3: "northern Argentina"}

INK, MUTED, SURFACE, GRID = "#1f2430", "#52514e", "#fbfaf7", "#9a9890"
CATEGORICAL = ["#2a78d6", "#eb6834", "#1f9e74"]  # blue, orange, green: one per group
BLUE = CATEGORICAL[0]
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


# ---------------------------------------------------------------------------------------------
# Reading one run
# ---------------------------------------------------------------------------------------------


def run_seconds(folder):
    """NAME's own total run time (s), from the timer summary in inputTiming.txt."""
    text = (folder / "inputTiming.txt").read_text()
    match = re.search(r"Timer TotalTime summary:\s*\n\s*Total time\s+([\d.]+)", text)
    if match is None:
        raise ValueError(f"No total time in {folder / 'inputTiming.txt'}")
    return float(match.group(1))


def output_megabytes(folder):
    """Size of the run's output folder (MB) and of its fields files alone."""
    files = list((folder / "output").iterdir())
    fields = [f for f in files if f.name.startswith("Fields_")]
    return (sum(f.stat().st_size for f in files) / 1e6,
            sum(f.stat().st_size for f in fields) / 1e6)  # fmt: skip


def run_cost(folder):
    """Run time, peak memory and output size of one run."""
    total_mb, fields_mb = output_megabytes(folder)
    return {
        "run_minutes": run_seconds(folder) / 60,
        "peak_memory_gb": int((folder / "peak_memory_kb.txt").read_text()) / 1e6,
        "output_mb": total_mb,
        "fields_mb": fields_mb,
    }


def group_weights(plan, obs):
    """{group number: (source names by level 1..16, column weights by level)}.

    Groups are numbered 1, 2, 3 in order of their IDs, as the per-group grids
    (prepare_unbounded_test.variant_grids). The weights are those of the sounding whose ID is the
    group ID (the group's earliest sounding).
    """
    sounding_ids = [str(s) for s in obs["sounding_id"].values]
    weights = (
        obs["pressure_weights"].values[:, :N_LEVELS]
        * obs["xch4_averaging_kernel"].values[:, :N_LEVELS]
    )
    out = {}
    for number, (group, rows) in enumerate(plan.groupby("release_group", sort=True), 1):
        rows = rows.sort_values("level_number")
        if list(rows["level_number"]) != list(range(1, N_LEVELS + 1)):
            raise ValueError(f"Group {group} does not have levels 1..{N_LEVELS}")
        out[number] = (list(rows["source_name"]), weights[sounding_ids.index(group)])
    return out


def column_footprint(fields, sources, weights):
    """Column footprint (lat, lon) from a fields DataArray: Σ_level w × footprint of the level."""
    area = grid_cell_area_m2(fields["lat"].values, fields["lon"].values)
    per_level = footprint_from_mixing_ratio(fields.sel(source=sources).values, area,
                                            RELEASE_SECONDS)  # fmt: skip
    total = np.tensordot(weights, per_level, axes=1)
    return xr.DataArray(total, dims=("lat", "lon"),
                        coords={"lat": fields["lat"], "lon": fields["lon"]})  # fmt: skip


def per_area(footprint):
    """Footprint per unit area (per km²): comparable between grids of different cell size."""
    area_km2 = grid_cell_area_m2(footprint["lat"].values, footprint["lon"].values) / 1e6
    return footprint / area_km2


def read_end_positions(folder, plan, weights_by_group):
    """Fate records (one per particle, X20) with group number and column weight per particle.

    Each particle carries weight w_level / (particles in its release), so a group's weights add
    up to its column weight, as the boundary sensitivities in processing.
    """
    particles, _ = read_name_fate_file(folder / "output" / "Particles_C1.txt")
    particles = drop_repeated_fate_records(particles, max_gap_seconds=SYNC_SECONDS)
    number_of = {}
    weight_of = {}
    for group, (sources, weights) in weights_by_group.items():
        for level, source in enumerate(sources):
            source_number = int(source.removeprefix("SourceID_"))
            number_of[source_number] = group
            weight_of[source_number] = weights[level]
    particles["group"] = particles["source_number"].map(number_of)
    counts = particles.groupby("source_number").size()
    particles["weight"] = particles["source_number"].map(weight_of) / particles[
        "source_number"
    ].map(counts)
    return particles


# ---------------------------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------------------------


def bar_panel(ax, labels, values, unit, fmt):
    """Horizontal bars, values labelled at the bar ends; recessive axes."""
    y = np.arange(len(labels))[::-1]
    ax.barh(y, values, height=0.6, color=BLUE)
    for yi, v in zip(y, values, strict=True):
        ax.text(v, yi, "  " + fmt.format(v), va="center", fontsize=16, color=INK)
    ax.set_yticks(y, labels)
    ax.set_xlabel(unit)
    ax.set_xlim(0, max(values) * 1.3)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


SMOKE_LABELS = {
    "bounded": "Today: bounded domain,\none grid",
    "per_group": "Unbounded, 0.1° grid\nper release group",
    "shared_region": "Unbounded, one 0.1° grid\nover the region",
    "shared_global": "Unbounded, one 0.1° grid\nover the globe",
}


def figure_smoke(out, smoke):
    names = list(SMOKE_LABELS)
    fig, axes = plt.subplots(1, 2, figsize=(18, 6.2), sharey=True)
    bar_panel(
        axes[0],
        [SMOKE_LABELS[n] for n in names],
        [smoke[n]["peak_memory_gb"] for n in names],
        "Peak memory (GB)",
        "{:.1f} GB",
    )
    bar_panel(axes[1], [SMOKE_LABELS[n] for n in names],
              [smoke[n]["run_minutes"] for n in names], "Run time, 1 h of travel (min)",
              "{:.1f} min")  # fmt: skip
    fig.tight_layout()
    fig.savefig(out / "fu1_smoke_cost.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


FULL_LABELS = {
    "bounded": "Today: bounded, 720 h",
    "per_group_240h": "Unbounded, 240 h",
    "per_group": "Unbounded, 720 h",
}


def figure_full(out, full):
    names = list(FULL_LABELS)
    fig, axes = plt.subplots(1, 2, figsize=(18, 4.8), sharey=True)
    bar_panel(axes[0], [FULL_LABELS[n] for n in names], [full[n]["run_minutes"] for n in names],
              "Run time, 48 releases, 4 threads (min)", "{:.0f} min")  # fmt: skip
    bar_panel(
        axes[1],
        [FULL_LABELS[n] for n in names],
        [full[n]["peak_memory_gb"] for n in names],
        "Peak memory (GB)",
        "{:.1f} GB",
    )
    fig.tight_layout()
    fig.savefig(out / "fu2_full_cost.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def add_domain_box(ax, **kwargs):
    d = DOMAIN
    ax.plot([d["x_min"], d["x_max"], d["x_max"], d["x_min"], d["x_min"]],
            [d["y_min"], d["y_min"], d["y_max"], d["y_max"], d["y_min"]],
            transform=ccrs.PlateCarree(), **kwargs)  # fmt: skip


def figure_end_positions(out, ends):
    """Where particles stop: today's exits vs unbounded end positions, coloured by group."""
    titles = {
        "bounded": "Today: where particles leave the domain",
        "per_group_240h": "Unbounded: positions after 240 h",
        "per_group": "Unbounded: positions after 720 h",
    }
    fig = plt.figure(figsize=(27, 5.4))
    for i, (name, title) in enumerate(titles.items()):
        ax = fig.add_subplot(1, 3, i + 1, projection=ccrs.Robinson(central_longitude=-40))
        ax.set_global()
        ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=0.6, edgecolor=MUTED)
        add_domain_box(ax, color=INK, linewidth=1.5, linestyle="--")
        p = ends[name]
        # Draw a random subsample, mixed across groups so no group hides another.
        sample = p.sample(n=min(len(p), 20000), random_state=1)
        for group, colour in zip(GROUP_NAMES, CATEGORICAL, strict=True):
            sel = sample[sample["group"] == group]
            ax.scatter(sel["lon"], sel["lat"], s=2, color=colour, alpha=0.5, linewidth=0,
                       transform=ccrs.PlateCarree(), rasterized=True)  # fmt: skip
        ax.set_title(title, loc="left", fontsize=20)
    handles = [plt.Line2D([], [], marker="o", linestyle="", color=c, markersize=10)
               for c in CATEGORICAL]  # fmt: skip
    fig.legend(handles, [f"Released: {n}" for n in GROUP_NAMES.values()],
               loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.02))  # fmt: skip
    fig.savefig(out / "fu3_end_positions.png", dpi=100, bbox_inches="tight")
    plt.close(fig)


def figure_nested(out, footprints, group):
    """Column footprint of one group on the nested grids, per unit area, log colour scale."""
    coarse = per_area(footprints["coarse"][group])
    medium = per_area(footprints["medium"][group])
    fine = per_area(footprints["fine"][group])
    vmax = float(fine.max())
    norm = mcolors.LogNorm(vmax * 1e-5, vmax)

    fig = plt.figure(figsize=(22, 8.5))
    layout = fig.add_gridspec(1, 2, width_ratios=[1.6, 1])
    ax = fig.add_subplot(layout[0], projection=ccrs.PlateCarree())
    ax.set_extent([-150, 40, -75, 50], ccrs.PlateCarree())
    for field, edge in ((coarse, None), (medium, INK), (fine, ORANGE)):
        mesh = ax.pcolormesh(
            field["lon"], field["lat"], field.where(field > 0), norm=norm, cmap=SEQUENTIAL,
            transform=ccrs.PlateCarree(), shading="nearest",
        )  # fmt: skip
        if edge is not None:
            lon, lat = field["lon"].values, field["lat"].values
            dx, dy = lon[1] - lon[0], lat[1] - lat[0]
            x0, x1, y0, y1 = lon[0] - dx / 2, lon[-1] + dx / 2, lat[0] - dy / 2, lat[-1] + dy / 2
            ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], color=edge, linewidth=1.5,
                    transform=ccrs.PlateCarree())  # fmt: skip
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.6, edgecolor=MUTED)
    ax.set_title("2.5° globe · 0.5° region (black) · 0.1° near the release (orange)",
                 loc="left", fontsize=17)  # fmt: skip

    ax2 = fig.add_subplot(layout[1], projection=ccrs.PlateCarree())
    lon, lat = fine["lon"].values, fine["lat"].values
    ax2.set_extent([lon[0], lon[-1], lat[0], lat[-1]], ccrs.PlateCarree())
    ax2.pcolormesh(lon, lat, fine.where(fine > 0), norm=norm, cmap=SEQUENTIAL,
                   transform=ccrs.PlateCarree(), shading="nearest")  # fmt: skip
    ax2.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.8, edgecolor=MUTED)
    ax2.add_feature(cfeature.RIVERS.with_scale("10m"), linewidth=0.5, edgecolor="#9ec5f4")
    ax2.set_title("Zoom: the 0.1° grid, ±5°", loc="left", fontsize=17)
    colorbar = fig.colorbar(mesh, ax=[ax, ax2], shrink=0.7, pad=0.02, extend="min")
    colorbar.set_label("Column footprint per area,\n(mol/mol) / (mol m⁻² s⁻¹) per km²")
    fig.savefig(out / "fu4_nested_footprint.png", dpi=100, bbox_inches="tight")
    plt.close(fig)


ORANGE = CATEGORICAL[1]


def inside_region(field):
    lon, lat = field["lon"], field["lat"]
    return ((lon > REGION["x_min"]) & (lon < REGION["x_max"]) & (lat > REGION["y_min"])
            & (lat < REGION["y_max"]))  # fmt: skip


def region_numbers(footprints_by_run, bounded):
    """Per group and run: footprint totals inside the region (medium grid) and outside it.

    Outside = the coarse global total minus the coarse cells inside the region. The test's region
    edges are not on 2.5° lines, so coarse cells along the edge are counted as inside or outside
    by their centre (a small error; a real tile scheme aligns every level).
    """
    rows = []
    for group, name in GROUP_NAMES.items():
        rows.append({"group": group, "name": name, "run": "bounded",
                     "inside": float(bounded[group].sum()), "outside": 0.0})  # fmt: skip
        for run, footprints in footprints_by_run.items():
            coarse = footprints["coarse"][group]
            coarse_inside = float(coarse.where(inside_region(coarse), 0).sum())
            rows.append(
                {
                    "group": group,
                    "name": name,
                    "run": run,
                    "inside": float(footprints["medium"][group].sum()),
                    "outside": float(coarse.sum()) - coarse_inside,
                }
            )
    return pd.DataFrame(rows)


REGION_RUNS = {
    "bounded": "Today (stops at the edge)",
    "per_group_240h": "Unbounded, 240 h",
    "per_group": "Unbounded, 720 h",
}


def figure_region(out, table):
    """Each group's column footprint total, inside and outside the region, by run."""
    fig, ax = plt.subplots(figsize=(16, 8))
    labels, y_positions = [], []
    y = 0.0
    for group, name in reversed(list(GROUP_NAMES.items())):
        for run in reversed(list(REGION_RUNS)):
            row = table[(table["group"] == group) & (table["run"] == run)].iloc[0]
            ax.barh(y, row["inside"], height=0.75, color=MUTED if run == "bounded" else BLUE)
            if row["outside"] > 0:
                ax.barh(y, row["outside"], left=row["inside"], height=0.75, color=ORANGE)
                share = row["outside"] / (row["inside"] + row["outside"])
                ax.text(row["inside"] + row["outside"], y, f"  {share:.0%} outside",
                        va="center", fontsize=15)  # fmt: skip
            labels.append(REGION_RUNS[run])
            y_positions.append(y)
            y += 1
        ax.text(-0.02, y - 0.35, name, transform=ax.get_yaxis_transform(), ha="right",
                fontsize=18, fontweight="bold")  # fmt: skip
        y += 0.9
    ax.set_yticks(y_positions, labels, fontsize=15)
    ax.set_xlabel("Total column footprint, (mol/mol) / (mol m⁻² s⁻¹)")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (MUTED, BLUE, ORANGE)]
    ax.legend(handles, ["Today: inside the domain", "Unbounded: inside the region",
                        "Unbounded: outside the region"],
              frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=3,
              fontsize=15)  # fmt: skip
    fig.tight_layout()
    fig.savefig(out / "fu5_region_vs_outside.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


DISTANCE_BANDS_KM = (0, 300, 1000, 3000, np.inf)
DISTANCE_LABELS = ("< 300 km", "300–1,000 km", "1,000–3,000 km", "> 3,000 km")


def distance_km(lat, lon, lat0, lon0):
    """Great-circle distance (km) of every grid cell centre (lat, lon) from (lat0, lon0)."""
    la, lo = np.meshgrid(np.radians(lat), np.radians(lon), indexing="ij")
    dlat, dlon = la - np.radians(lat0), lo - np.radians(lon0)
    a = np.sin(dlat / 2) ** 2 + np.cos(la) * np.cos(np.radians(lat0)) * np.sin(dlon / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def distance_numbers(footprints_by_run, bounded, plan):
    """Per group and run: column footprint total in each distance band from the release.

    Uses today's grid for the bounded run and the 0.5° regional grid for the unbounded runs, so
    only cells inside the region count. The two grids differ, so bands of a few cells (< 100 km)
    would differ by binning alone; the smallest band is 300 km.
    """
    rows = []
    for number, (_group, releases) in enumerate(plan.groupby("release_group", sort=True), 1):
        lat0, lon0 = releases["lat"].iloc[0], releases["lon"].iloc[0]
        runs = {"bounded": bounded[number]}
        runs.update({run: fps["medium"][number] for run, fps in footprints_by_run.items()})
        for run, fp in runs.items():
            d = distance_km(fp["lat"].values, fp["lon"].values, lat0, lon0)
            for low, high, label in zip(DISTANCE_BANDS_KM[:-1], DISTANCE_BANDS_KM[1:],
                                        DISTANCE_LABELS, strict=True):  # fmt: skip
                total = float(fp.values[(d >= low) & (d < high)].sum())
                rows.append({"group": number, "name": GROUP_NAMES[number], "run": run,
                             "band": label, "total": total})  # fmt: skip
    return pd.DataFrame(rows)


def figure_distance(out, table):
    """Column footprint by distance band, one panel per group, today / 240 h / 720 h."""
    colours = {"bounded": MUTED, "per_group_240h": BLUE, "per_group": ORANGE}
    fig, axes = plt.subplots(1, 3, figsize=(24, 6.5), sharey=False)
    x = np.arange(len(DISTANCE_LABELS))
    width = 0.27
    for ax, (group, name) in zip(axes, GROUP_NAMES.items(), strict=True):
        for i, run in enumerate(REGION_RUNS):
            sel = table[(table["group"] == group) & (table["run"] == run)].set_index("band")
            ax.bar(x + (i - 1) * width, sel.loc[list(DISTANCE_LABELS), "total"], width=width,
                   color=colours[run], label=REGION_RUNS[run])  # fmt: skip
        ax.set_xticks(x, [label.replace(" km", "\nkm") for label in DISTANCE_LABELS], fontsize=15)
        ax.set_title(name, loc="left", fontsize=20)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("Column footprint inside the region,\n(mol/mol) / (mol m⁻² s⁻¹)")
    axes[0].legend(frameon=False, fontsize=15, loc="upper left")
    fig.tight_layout()
    fig.savefig(out / "fu6_by_distance.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------------


def read_grid_fields(folder, grid_name):
    return read_fields(next((folder / "output").glob(f"Fields_{grid_name}_C1_T1_*.txt")))


def unbounded_footprints(run, weights):
    """{"fine"|"medium"|"coarse": {group: column footprint}} of one unbounded run."""
    medium = read_grid_fields(run, "HGridMedium")
    coarse = read_grid_fields(run, "HGridCoarse")
    return {
        "medium": {g: column_footprint(medium, *weights[g]) for g in weights},
        "coarse": {g: column_footprint(coarse, *weights[g]) for g in weights},
        "fine": {
            g: column_footprint(read_grid_fields(run, f"HGridFineG{g}"), *weights[g])
            for g in weights
        },
    }


def main(smoke_folder, full_folder, obs_file, out):
    smoke_folder, full_folder, out = Path(smoke_folder), Path(full_folder), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    obs = xr.open_dataset(obs_file).load()

    smoke = {name: run_cost(smoke_folder / name) for name in SMOKE_LABELS}
    full = {name: run_cost(full_folder / name) for name in FULL_LABELS}

    plan = pd.read_csv(full_folder / "per_group" / "release_plan.csv",
                       dtype={"release_group": str})  # fmt: skip
    weights = group_weights(plan, obs)

    footprints_by_run = {
        run: unbounded_footprints(full_folder / run, weights)
        for run in ("per_group_240h", "per_group")
    }
    footprints = footprints_by_run["per_group"]
    bounded_fields = read_grid_fields(full_folder / "bounded", "HGrid1")
    bounded = {g: column_footprint(bounded_fields, *weights[g]) for g in weights}

    ends = {name: read_end_positions(full_folder / name, plan, weights) for name in FULL_LABELS}
    end_summary = {}
    for name, p in ends.items():
        end_summary[name] = {
            "records": len(p),
            "fate_codes": {str(k): int(v) for k, v in p["fate_code"].value_counts().items()},
            "mean_travel_hours": float(p["age_hours"].mean()),
            "median_height_m_agl": float(p["z"].median()),
            "share_south_of_60S": float(p.loc[p["lat"] < -60, "weight"].sum() / p["weight"].sum()),
            "share_north_of_equator": float(
                p.loc[p["lat"] > 0, "weight"].sum() / p["weight"].sum()
            ),
        }

    table = region_numbers(footprints_by_run, bounded)
    by_distance = distance_numbers(footprints_by_run, bounded, plan)
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                            text=True, check=False).stdout.strip()  # fmt: skip
    summary = {"smoke": smoke, "full": full, "end_positions": end_summary,
               "region": table.to_dict(orient="records"),
               "by_distance": by_distance.to_dict(orient="records"), "commit": commit}  # fmt: skip
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    figure_smoke(out, smoke)
    figure_full(out, full)
    figure_end_positions(out, ends)
    figure_nested(out, footprints, group=2)
    figure_region(out, table)
    figure_distance(out, by_distance)


if __name__ == "__main__":
    main(*sys.argv[1:5])
