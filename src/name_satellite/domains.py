"""Domains and satellite sites: where NAME runs and which soundings it uses.

Where this sits (see docs/workflow.md)
--------------------------------------
Configuration used by steps 1 (which soundings: ``Site``) and 3 (the NAME grids: ``Domain``).

A ``Domain`` is the NAME output grid (where footprints are reported) and the dispersion domain
(where particles are followed; they stop when they leave it). A ``Site`` is a satellite selection
for one domain: the lat/lon box of soundings to use and, optionally, legacy binning.

Values are those of the legacy 2013 GOSAT-BRAZIL runs unless noted; each cites where it came from.
Add new domains and sites here rather than in scripts.
"""

from dataclasses import dataclass

from name_satellite.name_input import DispersionDomain, OutputGrid


@dataclass(frozen=True)
class Domain:
    """A NAME domain: output grid and dispersion domain (degrees; domain top in m agl)."""

    name: str
    grid: OutputGrid
    dispersion: DispersionDomain


@dataclass(frozen=True)
class Site:
    """A satellite selection within a domain.

    Attributes:
        name: site name as used in file names (legacy "GOSAT-BRAZIL").
        domain: the ``Domain`` its NAME runs use.
        lat_bounds, lon_bounds: soundings with lower <= value < upper are used (degrees).
        bin_size: (lat, lon) degrees for ``observations.legacy_bin``, or None for no binning.
        footprint_km: diameter of a sounding's ground footprint, for the release box.
    """

    name: str
    domain: Domain
    lat_bounds: tuple
    lon_bounds: tuple
    bin_size: tuple | None = None
    footprint_km: float = 10.5


# South America: output grid of the 2019 mixing-ratio run (legacy input file
# BackRun_SOUTHAMERICA_GOSAT-BRAZIL_column_20130101.txt, "Horizontal Grids" HGrid1: 190 x 357 cells
# of 0.352 x 0.234 degrees, first cell centre -91.329, -60.981) and its dispersion domain (same
# file, "Domains", copied verbatim): 4.0-5.4 degrees beyond the outer cell centres (W 5.0, E 5.35,
# S 4.0, N 5.23), top 19 km agl.
SOUTHAMERICA = Domain(
    name="SOUTHAMERICA",
    grid=OutputGrid(nx=190, ny=357, dx=0.352, dy=0.234, x_min=-91.329, y_min=-60.981),
    dispersion=DispersionDomain(x_min=-96.329, x_max=-19.449, y_min=-64.981, y_max=27.557),
)

# GOSAT over Brazil: lat/lon box from the 2013 legacy observation file's "bins" attribute.
# Legacy binned the soundings to the NAME grid (0.234 x 0.352); new runs use single soundings
# (bin_size None) until the group decides otherwise (decisions.md open question 6).
GOSAT_BRAZIL = Site(
    name="GOSAT-BRAZIL",
    domain=SOUTHAMERICA,
    lat_bounds=(-35.753, 7.251),
    lon_bounds=(-75.984, -32.793),
    bin_size=None,
)

SITES = {site.name: site for site in (GOSAT_BRAZIL,)}
