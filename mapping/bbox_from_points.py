"""
vibecoded with Claude

Generate a bounding box around two lat/lon points, with UTM support.

Takes two lat/lon points (e.g. an event and a station, or any two corners
of interest), buffers outward by a user-controlled distance, and reports
the four corners of the resulting bounding box in both lon/lat and UTM.

Buffering is done in UTM (meters) so the requested distance is accurate,
then the result is converted back to lon/lat.

Example
-------
    python bbox_from_points.py --p1 61.2 -149.9 --p2 61.5 -150.3 --buffer 25

    python bbox_from_points.py --p1 61.2 -149.9 --p2 61.5 -150.3 \\
        --buffer-north 10 --buffer-south 30 --buffer-east 15 --buffer-west 15
"""
import argparse
from dataclasses import dataclass
from typing import Optional, Tuple

from pyproj import Transformer


@dataclass
class BoundingBox:
    min_lon: float
    max_lon: float
    min_lat: float
    max_lat: float
    zone: int
    hemisphere: str
    min_easting: float
    max_easting: float
    min_northing: float
    max_northing: float

    def corners_lonlat(self) -> dict:
        return {
            "NW": (self.min_lon, self.max_lat),
            "NE": (self.max_lon, self.max_lat),
            "SE": (self.max_lon, self.min_lat),
            "SW": (self.min_lon, self.min_lat),
        }

    def corners_utm(self) -> dict:
        return {
            "NW": (self.min_easting, self.max_northing),
            "NE": (self.max_easting, self.max_northing),
            "SE": (self.max_easting, self.min_northing),
            "SW": (self.min_easting, self.min_northing),
        }


def utm_zone_from_lonlat(lon: float, lat: float) -> Tuple[int, str]:
    """Return the (zone, hemisphere) a lon/lat point falls into."""
    zone = int((lon + 180) / 6) + 1
    hemisphere = "N" if lat >= 0 else "S"
    return zone, hemisphere


def utm_crs(zone: int, hemisphere: str) -> str:
    return f"EPSG:{32600 + zone if hemisphere == 'N' else 32700 + zone}"


def compute_bounding_box(
    p1: Tuple[float, float],
    p2: Tuple[float, float],
    buffer_km: float = 0.0,
    buffer_north_km: Optional[float] = None,
    buffer_south_km: Optional[float] = None,
    buffer_east_km: Optional[float] = None,
    buffer_west_km: Optional[float] = None,
    zone: Optional[int] = None,
    hemisphere: Optional[str] = None,
) -> BoundingBox:
    """Compute a bounding box around two (lat, lon) points.

    Parameters
    ----------
    p1, p2
        (lat, lon) points to bound.
    buffer_km
        Default distance (km) added to every edge. Overridden per-edge by
        the buffer_<direction>_km arguments below.
    buffer_north_km, buffer_south_km, buffer_east_km, buffer_west_km
        Optional per-edge distances (km), overriding buffer_km for that edge.
    zone, hemisphere
        UTM zone/hemisphere to project into. Auto-detected from the
        midpoint of p1/p2 if not given.
    """
    lat1, lon1 = p1
    lat2, lon2 = p2

    if zone is None or hemisphere is None:
        mid_lon = (lon1 + lon2) / 2
        mid_lat = (lat1 + lat2) / 2
        zone, hemisphere = utm_zone_from_lonlat(mid_lon, mid_lat)

    crs = utm_crs(zone, hemisphere)
    to_utm = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    to_lonlat = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)

    e1, n1 = to_utm.transform(lon1, lat1)
    e2, n2 = to_utm.transform(lon2, lat2)

    north_m = (buffer_north_km if buffer_north_km is not None else buffer_km) * 1000
    south_m = (buffer_south_km if buffer_south_km is not None else buffer_km) * 1000
    east_m = (buffer_east_km if buffer_east_km is not None else buffer_km) * 1000
    west_m = (buffer_west_km if buffer_west_km is not None else buffer_km) * 1000

    min_easting = min(e1, e2) - west_m
    max_easting = max(e1, e2) + east_m
    min_northing = min(n1, n2) - south_m
    max_northing = max(n1, n2) + north_m

    min_lon, min_lat = to_lonlat.transform(min_easting, min_northing)
    max_lon, max_lat = to_lonlat.transform(max_easting, max_northing)

    return BoundingBox(
        min_lon=min_lon,
        max_lon=max_lon,
        min_lat=min_lat,
        max_lat=max_lat,
        zone=zone,
        hemisphere=hemisphere,
        min_easting=min_easting,
        max_easting=max_easting,
        min_northing=min_northing,
        max_northing=max_northing,
    )


def print_bounding_box(bbox: BoundingBox) -> None:
    print(f"UTM zone: {bbox.zone}{bbox.hemisphere}")
    print()
    print("Corners (lon, lat):")
    for name, (lon, lat) in bbox.corners_lonlat().items():
        print(f"  {name}: {lon:.6f}, {lat:.6f}")
    print()
    print("Corners (easting, northing) [m]:")
    for name, (e, n) in bbox.corners_utm().items():
        print(f"  {name}: {e:.1f}, {n:.1f}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p1", nargs=2, type=float, required=True, metavar=("LAT", "LON"))
    parser.add_argument("--p2", nargs=2, type=float, required=True, metavar=("LAT", "LON"))
    parser.add_argument("--buffer", type=float, default=0.0, help="Default buffer distance (km) applied to all edges.")
    parser.add_argument("--buffer-north", type=float, default=None, help="Override buffer (km) on the north edge.")
    parser.add_argument("--buffer-south", type=float, default=None, help="Override buffer (km) on the south edge.")
    parser.add_argument("--buffer-east", type=float, default=None, help="Override buffer (km) on the east edge.")
    parser.add_argument("--buffer-west", type=float, default=None, help="Override buffer (km) on the west edge.")
    parser.add_argument("--zone", type=int, default=None, help="UTM zone to use (auto-detected if omitted).")
    parser.add_argument("--hemisphere", choices=["N", "S"], default=None, help="UTM hemisphere (auto-detected if omitted).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    bbox = compute_bounding_box(
        p1=tuple(args.p1),
        p2=tuple(args.p2),
        buffer_km=args.buffer,
        buffer_north_km=args.buffer_north,
        buffer_south_km=args.buffer_south,
        buffer_east_km=args.buffer_east,
        buffer_west_km=args.buffer_west,
        zone=args.zone,
        hemisphere=args.hemisphere,
    )
    print_bounding_box(bbox)


if __name__ == "__main__":
    main()
