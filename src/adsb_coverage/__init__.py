"""Theoretical reception coverage of a fixed ADS-B / ADS-L ground receiver (line of sight, terrain, k-factor)."""
from .core import (  # noqa: F401
    GeoTiffMosaic, SwissAlti3D, Terrain, compute_coverage, continuous_range, copernicus_urls,
    decimate_xyz_zip, destination, main, max_visible_range, visible_runs, wgs84_to_lv95, write_kml,
)

__version__ = "0.1.0"
