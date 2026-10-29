"""Couverture de réception théorique d'un récepteur ADS-B / ADS-L fixe (LOS, relief, k-factor)."""
from .core import (  # noqa: F401
    SwissAlti3D, GeoTiffMosaic, Terrain, compute_coverage, continuous_range,
    max_visible_range, coverage_ranges, visible_runs, write_kml, wgs84_to_lv95, destination, main,
)

__version__ = "0.1.0"
