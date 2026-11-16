"""The GLIDE reference-run specification, as numbers: receptors, domains, blocks and NAME settings.

Where this sits
---------------
Configuration only (no files read or written). The source is GLIDE's
``docs/reference-runs.md`` (repository github.com/mrghg/glide); every value below cites its
section ("§"). The version used is recorded in ``SPEC_VERSION``, and goes into the run record.
``plan.py`` turns the receptors and blocks into a release plan; ``prepare.py`` writes the NAME
input with ``reference_settings``.

These runs are separate from the satellite pipeline: they share its NAME input writer
(``name_input``), job scripts (``lotus``) and met handling (``met``), but none of its sounding or
column code.

Rules of the spec that shape this file (§ preamble): use its numbers exactly; keep every setting
the same in all runs of one variant; settings that NAME and FLEXPART both have take the same
value, even where that differs from the group's usual NAME setup.
"""

from dataclasses import dataclass, replace

import pandas as pd

from name_satellite.domains import Domain
from name_satellite.met import MET_VERSIONS, NAME_DIR
from name_satellite.name_input import DispersionDomain, OutputGrid, RunSettings

# GLIDE repository commit of docs/reference-runs.md these numbers come from (branch
# docs/reference-run-settings, PR #21).
SPEC_VERSION = "glide e934448 docs/reference-runs.md"

# ---------------------------------------------------------------------------------------------
# Domains and output grids (§2)
# ---------------------------------------------------------------------------------------------
# The particle domain is the box; particles are removed when they leave it. The output grid's
# cells are 0.25 degrees, centred on ERA5 grid points, and stop half a cell inside the box.
# OutputGrid.x_min / y_min are the FIRST CELL CENTRES (NAME's X Min, Y Min of a grid).
# Top: 20,000 m above ground (§2: no particle may be removed below 20 km above the ground).
DOMAIN_TOP_M_AGL = 20000.0

EUR = Domain(
    name="EUR",
    grid=OutputGrid(nx=279, ny=159, dx=0.25, dy=0.25, x_min=-29.75, y_min=32.25),
    dispersion=DispersionDomain(
        x_min=-30.0, x_max=40.0, y_min=32.0, y_max=72.0, z_max_m_agl=DOMAIN_TOP_M_AGL
    ),
)
NAM = Domain(
    name="NAM",
    grid=OutputGrid(nx=259, ny=143, dx=0.25, dy=0.25, x_min=-119.75, y_min=22.25),
    dispersion=DispersionDomain(
        x_min=-120.0, x_max=-55.0, y_min=22.0, y_max=58.0, z_max_m_agl=DOMAIN_TOP_M_AGL
    ),
)
DOMAINS = {domain.name: domain for domain in (EUR, NAM)}


# ---------------------------------------------------------------------------------------------
# Receptors (§3)
# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Receptor:
    """One receptor. Heights are metres above the MODEL's ground (§3), never above sea level.

    ``height_m_agl`` is None for the column release at Lamont (§5).
    """

    code: str
    site: str
    lat: float  # degrees north
    lon: float  # degrees east (west is negative)
    height_m_agl: float | None
    domain: str  # "EUR" or "NAM"

    @property
    def is_column(self):
        return self.height_m_agl is None


RECEPTORS = {
    r.code: r
    for r in (
        Receptor("MHD", "Mace Head", 53.3267, -9.9046, 10.0, "EUR"),
        Receptor("JFJ", "Jungfraujoch", 46.5475, 7.9851, 1000.0, "EUR"),
        Receptor("CBW-20", "Cabauw, low", 51.9703, 4.9264, 20.0, "EUR"),
        Receptor("CBW-200", "Cabauw, high", 51.9703, 4.9264, 200.0, "EUR"),
        Receptor("TAC", "Tacolneston", 52.5176, 1.1387, 185.0, "EUR"),
        Receptor("LEF", "Park Falls", 45.9451, -90.2732, 30.0, "NAM"),
        Receptor("SGP", "Lamont", 36.6040, -97.4860, None, "NAM"),
    )
}


# ---------------------------------------------------------------------------------------------
# Blocks: release times and the span the run must cover (§4)
# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Block:
    """Two weeks of hourly releases. Times are UTC; a release is named by its start time."""

    name: str  # e.g. "EUR-Jan"
    domain: str
    first_release: pd.Timestamp
    last_release: pd.Timestamp  # start of the last release
    cover_start: pd.Timestamp  # "Run must cover" (§4)
    cover_end: pd.Timestamp

    @property
    def month(self):
        """YYYYMM, as used in the delivered file names (§8)."""
        return self.first_release.strftime("%Y%m")

    def release_starts(self):
        """Start times of every hourly release (336 per block)."""
        return pd.date_range(self.first_release, self.last_release, freq="1h")


def _block(name, domain, month):
    jan = month == "Jan"
    year_month = "2019-01" if jan else "2019-07"
    return Block(
        name=name,
        domain=domain,
        first_release=pd.Timestamp(f"{year_month}-01 00:00"),
        last_release=pd.Timestamp(f"{year_month}-14 23:00"),
        cover_start=pd.Timestamp("2018-12-22 00:00" if jan else "2019-06-21 00:00"),
        cover_end=pd.Timestamp(f"{year_month}-15 00:00"),
    )


BLOCKS = {
    b.name: b
    for b in (
        _block("EUR-Jan", "EUR", "Jan"),
        _block("EUR-Jul", "EUR", "Jul"),
        _block("NAM-Jan", "NAM", "Jan"),
        _block("NAM-Jul", "NAM", "Jul"),
    )
}

# ---------------------------------------------------------------------------------------------
# Releases (§1, §4, §5)
# ---------------------------------------------------------------------------------------------
RELEASE_HOURS = 1  # each release lasts one hour, particles evenly through it (§4)
MAX_AGE_HOURS = 240  # each particle is followed for 240 h from its own release (§4)
# NAME takes a release rate, not a total (§4): 6/s = 21,600 per hour for a point,
# 28/s = 100,800 for the Lamont column. The number actually released is counted from the output.
POINT_PARTICLES_PER_SECOND = 6
COLUMN_PARTICLES_PER_SECOND = 28
COLUMN_TOP_HPA = 100.0  # the Lamont column runs from the model surface pressure to 100 hPa (§5)

# ---------------------------------------------------------------------------------------------
# Met (§1, §7): UM global analyses, Mk10 series, one series per block.
# ---------------------------------------------------------------------------------------------
MET_VERSION = next(v for v in MET_VERSIONS if v.name == "UMG_Mk10PT")
# NAME's met definition files (read, never copied), for the run record.
MET_VERSION_DEFINITIONS = NAME_DIR / "Resources" / "Defns"

# ---------------------------------------------------------------------------------------------
# NAME settings (§7, NAME table)
# ---------------------------------------------------------------------------------------------
SYNC_MINUTES = 5  # 300 s: NAME's main time step and footprint sampling interval (§7)
MIN_BL_DEPTH_M = 100.0  # boundary-layer depth limited to 100-4,500 m (§7)
MAX_BL_DEPTH_M = 4500.0
DEEP_CONVECTION = "New"  # NAME's mass-flux scheme (§7)
# Max # Particles: at least 6,000,000 for a block run of one point receptor (240 live releases ×
# 21,600 = 5.2 M), 25,000,000 for Lamont (§7, "The particle cap").
MAX_PARTICLES_POINT = 6_000_000
MAX_PARTICLES_COLUMN = 25_000_000
RANDOM_SEED_MAIN = "Fixed"  # `main` and the single-release runs (§7, §9)
RANDOM_SEED_REPEATS = "Random"  # seed1 ... seed5 (§9)


def reference_settings(run_name, output_folder, domain, met_folder, max_particles, threads=4):
    """NAME ``RunSettings`` for a `main` reference run (§7).

    Args:
        run_name: NAME's run name.
        output_folder: where NAME writes its output.
        domain: ``EUR`` or ``NAM`` (a ``Domain`` of this module).
        met_folder: where the met is restored (scratch).
        max_particles: NAME's cap on particles alive at once (MAX_PARTICLES_POINT or _COLUMN).
        threads: OpenMP threads (= Slurm CPUs).

    Returns:
        ``RunSettings``. Variants (`group`, `nomeso`, repeat seeds) are made from it with
        ``variant_settings``.
    """
    return RunSettings(
        run_name=run_name,
        output_folder=str(output_folder),
        grid=domain.grid,
        domain=domain.dispersion,
        met_version=MET_VERSION,
        met_folder=str(met_folder),
        max_age_hours=MAX_AGE_HOURS,
        sync_minutes=SYNC_MINUTES,
        max_particles=max_particles,
        threads=threads,
        random_seed=RANDOM_SEED_MAIN,
        min_bl_depth_m=MIN_BL_DEPTH_M,
        max_bl_depth_m=MAX_BL_DEPTH_M,
        deep_convection=DEEP_CONVECTION,
        mesoscale_motions=True,
        surface_layer_depth_m=40.0,  # footprint layer 0-40 m above the model's ground (§6)
        fate_z_coord="m agl",
        fate_met=False,
        integration="all_releases",
        release_met=False,  # sounding-based; receptor met is requested separately
    )


# The run variants of §9, each as the changes it makes to `main`'s settings.
VARIANTS = {
    "main": {},
    # Repeat-seed runs: NAME has no numeric seed, so each is a separate `Random` run (§9).
    **{f"seed{i}": {"random_seed": RANDOM_SEED_REPEATS} for i in range(1, 6)},
    # Unresolved mesoscale motions off; the seven standard cases (§9).
    "nomeso": {"mesoscale_motions": False},
    # The group's site template: deep convection off, 3-min sync, 40-4,000 m (§9).
    "group": {
        "deep_convection": "No",
        "sync_minutes": 3,
        "min_bl_depth_m": 40.0,
        "max_bl_depth_m": 4000.0,
    },
}


def variant_settings(settings, variant):
    """``settings`` with the changes of run variant ``variant`` (a key of ``VARIANTS``)."""
    if variant not in VARIANTS:
        raise ValueError(f"Unknown run variant {variant!r}; choose from {sorted(VARIANTS)}")
    return replace(settings, **VARIANTS[variant])
