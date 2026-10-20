"""Create a PyGMT map for a bounding box and add markers for stations/cities.

This script plots topography, contour elevations, coastlines, and allows
simple latitude/longitude markers with optional labels. Supports both lon/lat
and UTM coordinate systems.

* Vibecoded by GitHub Copilot
"""
import argparse
import contextlib
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import xml.etree.ElementTree as ET
from typing import List, Dict, Optional, Tuple

import pygmt
from pyproj import Transformer


def detect_coordinate_system(
    values: List[float] | Dict[str, float],
    utm_zone: Optional[int] = None,
) -> str:
    """Detect whether coordinates are in UTM or lon/lat format.

    Parameters
    ----------
    values
        Either a list of values [minx, maxx, miny, maxy] or a dict with keys
        for coordinates (lon/lat or easting/northing).
    utm_zone
        If provided, indicates UTM system. If None, detection is based on value ranges.

    Returns
    -------
    str
        Either "utm" or "lonlat".
    """
    if utm_zone is not None:
        return "utm"

    if isinstance(values, dict):
        if "easting" in values or "northing" in values:
            return "utm"
        if "lon" in values or "lat" in values:
            return "lonlat"

    if isinstance(values, (list, tuple)):
        # Check ranges: lon/lat are bounded, UTM values are typically larger
        # lon: -180 to 180, lat: -90 to 90
        # UTM easting: ~150,000 to ~850,000, northing: 0 to 10,000,000
        if len(values) >= 2:
            max_val = max(abs(float(v)) for v in values[:4])
            if max_val > 180:  # Likely UTM
                return "utm"
            else:  # Likely lon/lat
                return "lonlat"

    return "lonlat"  # Default


def utm_to_lonlat(
    easting: float, northing: float, zone: int, hemisphere: str = "N"
) -> Tuple[float, float]:
    """Convert UTM coordinates to lon/lat.

    Parameters
    ----------
    easting
        UTM easting coordinate.
    northing
        UTM northing coordinate.
    zone
        UTM zone number (1-60).
    hemisphere
        "N" for Northern hemisphere, "S" for Southern hemisphere.

    Returns
    -------
    tuple
        (lon, lat)
    """
    crs_utm = f"EPSG:{32600 + zone if hemisphere == 'N' else 32700 + zone}"
    transformer = Transformer.from_crs(crs_utm, "EPSG:4326", always_xy=True)
    lon, lat = transformer.transform(easting, northing)
    return lon, lat


def lonlat_to_utm(lon: float, lat: float) -> Tuple[float, float, int, str]:
    """Convert lon/lat to UTM coordinates.

    Parameters
    ----------
    lon
        Longitude.
    lat
        Latitude.

    Returns
    -------
    tuple
        (easting, northing, zone, hemisphere)
    """
    # Determine zone and hemisphere
    zone = int((lon + 180) / 6) + 1
    hemisphere = "N" if lat >= 0 else "S"

    crs_utm = f"EPSG:{32600 + zone if hemisphere == 'N' else 32700 + zone}"
    transformer = Transformer.from_crs("EPSG:4326", crs_utm, always_xy=True)
    easting, northing = transformer.transform(lon, lat)
    return easting, northing, zone, hemisphere


def convert_region_to_lonlat(
    region: List[float], utm_zone: Optional[int] = None, hemisphere: str = "N"
) -> List[float]:
    """Convert a region bounding box from UTM to lon/lat.

    Parameters
    ----------
    region
        [min_easting, max_easting, min_northing, max_northing] for UTM
        or [minlon, maxlon, minlat, maxlat] for lon/lat.
    utm_zone
        UTM zone number. If provided, assumes region is in UTM.
    hemisphere
        "N" for Northern, "S" for Southern hemisphere.

    Returns
    -------
    list
        [minlon, maxlon, minlat, maxlat]
    """
    if utm_zone is None:
        # Already in lon/lat
        return region

    min_easting, max_easting, min_northing, max_northing = region[:4]

    # Convert corners to determine actual lon/lat bounds
    corners_utm = [
        (min_easting, min_northing),
        (max_easting, min_northing),
        (min_easting, max_northing),
        (max_easting, max_northing),
    ]

    lons = []
    lats = []
    for easting, northing in corners_utm:
        lon, lat = utm_to_lonlat(easting, northing, utm_zone, hemisphere)
        lons.append(lon)
        lats.append(lat)

    return [min(lons), max(lons), min(lats), max(lats)]


def convert_markers_to_lonlat(
    markers: List[Dict[str, object]],
    utm_zone: Optional[int] = None,
    hemisphere: str = "N",
) -> List[Dict[str, object]]:
    """Convert marker coordinates from UTM to lon/lat if needed.

    Parameters
    ----------
    markers
        List of marker dictionaries. Each marker should have either
        (lon, lat) or (easting, northing) keys.
    utm_zone
        UTM zone number. If provided, assumes markers are in UTM.
    hemisphere
        "N" for Northern, "S" for Southern hemisphere.

    Returns
    -------
    list
        Markers with lon/lat coordinates.
    """
    converted_markers = []

    for marker in markers:
        marker_copy = marker.copy()

        if utm_zone is not None and "easting" in marker_copy and "northing" in marker_copy:
            # Convert from UTM to lon/lat
            easting = float(marker_copy.pop("easting"))
            northing = float(marker_copy.pop("northing"))
            lon, lat = utm_to_lonlat(easting, northing, utm_zone, hemisphere)
            marker_copy["lon"] = lon
            marker_copy["lat"] = lat

        converted_markers.append(marker_copy)

    return converted_markers


def get_region_area(region: List[float]) -> float:
    """Calculate the area of a region bounding box.

    Parameters
    ----------
    region
        [minx, maxx, miny, maxy]

    Returns
    -------
    float
        Area (width * height)
    """
    minx, maxx, miny, maxy = region[:4]
    width = maxx - minx
    height = maxy - miny
    return width * height


def find_largest_region(regions: List[List[float]]) -> Tuple[List[float], List[List[float]]]:
    """Find the largest region and return it along with the remaining smaller regions.

    Parameters
    ----------
    regions
        List of region bounding boxes, each [minx, maxx, miny, maxy]

    Returns
    -------
    tuple
        (largest_region, smaller_regions)
    """
    if len(regions) == 1:
        return regions[0], []
    
    areas = [get_region_area(r) for r in regions]
    max_idx = areas.index(max(areas))
    largest = regions[max_idx]
    smaller = [regions[i] for i in range(len(regions)) if i != max_idx]
    return largest, smaller


def nice_interval(span: float, target: int = 3) -> float:
    """Choose a round decimal-degree annotation interval for an axis.

    :type span: float
    :param span: axis extent in degrees
    :type target: int
    :param target: approximate number of annotations wanted
    :rtype: float
    :return: the largest interval from a 1-2-2.5-5 series giving at least
        `target` annotations, or the smallest series value if none does
    """
    candidates = [0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5,
                  1, 2, 2.5, 5, 10, 20, 30]
    for step in reversed(candidates):
        if span / step >= target:
            return step
    return candidates[0]


# Public Overpass API servers, tried in turn. These are shared, often
# overloaded services: when busy they answer HTTP 429/502/503/504 or simply
# never respond, so rotate through independent servers instead of waiting
# on one. (overpass.kumi.systems is an alias of overpass.private.coffee.)
OVERPASS_URLS = [
    "https://overpass-api.de/api/interpreter",
    # "https://maps.mail.ru/osm/tools/overpass/api/interpreter",
    "https://overpass.private.coffee/api/interpreter",
]

# Server-side limits sent with each query. Busy servers admit queries that
# ask for less time and memory more readily; ours are small bounding boxes
OVERPASS_SETTINGS = "[out:json][timeout:60][maxsize:134217728]"

# Seconds to wait for a response before treating a server as busy
OVERPASS_READ_TIMEOUT = 90


def overpass_query(query: str, rounds: int = 2,
                   wait: float = 10.0) -> List[Dict[str, object]]:
    """Run an Overpass API query, rotating through busy servers.

    Each round tries every server in OVERPASS_URLS once. A server that
    answers with a temporary error (HTTP 429, 502, 503, 504), fails to
    connect or does not respond within OVERPASS_READ_TIMEOUT seconds is
    skipped for the rest of that round. Between rounds the function waits
    `wait` seconds, doubling each time. Worst case it takes about
    rounds * len(OVERPASS_URLS) * OVERPASS_READ_TIMEOUT seconds.

    :type query: str
    :param query: Overpass QL query requesting JSON output
    :type rounds: int
    :param rounds: number of passes over all servers
    :type wait: float
    :param wait: seconds to wait after the first full round of failures
    :rtype: list of dict
    :return: the 'elements' list of the Overpass JSON response
    :raises RuntimeError: if every server fails in every round
    """
    data = urllib.parse.urlencode({"data": query}).encode()
    headers = {"User-Agent": "plot_bbox_map.py (pygmt station maps)"}
    errors = []
    delay = wait
    for round_ in range(rounds):
        for url in OVERPASS_URLS:
            host = urllib.parse.urlparse(url).netloc
            print(f"  querying {host} ...", flush=True)
            request = urllib.request.Request(url, data=data,
                                             headers=headers)
            try:
                with urllib.request.urlopen(
                        request, timeout=OVERPASS_READ_TIMEOUT) as resp:
                    return json.load(resp)["elements"]
            except urllib.error.HTTPError as e:
                if e.code not in (429, 502, 503, 504):
                    raise  # not a busy server, e.g. a malformed query
                msg = f"HTTP {e.code} {e.reason}"
            except (urllib.error.URLError, TimeoutError, OSError,
                    json.JSONDecodeError) as e:
                msg = str(getattr(e, "reason", e))
            print(f"  {host}: {msg}")
            errors.append(f"{host}: {msg}")
        if round_ < rounds - 1:
            print(f"  all servers busy; retrying in {delay:.0f} s")
            time.sleep(delay)
            delay *= 2
    raise RuntimeError(
        "all Overpass servers are busy or unreachable; try again later. "
        "Once a fetch succeeds it is cached, so this is only needed once "
        "per region. Errors:\n  " + "\n  ".join(errors))


# OpenStreetMap 'highway' classes fetched by default, from largest to smallest
DEFAULT_ROAD_TYPES = ["motorway", "trunk", "primary", "secondary", "tertiary"]

# OpenStreetMap 'waterway' classes fetched by default
DEFAULT_RIVER_TYPES = ["river"]


def fetch_osm_lines(
    region: List[float],
    tag: str,
    types: List[str],
    cache: Optional[str] = None,
) -> List[Dict[str, object]]:
    """Fetch OpenStreetMap line features of given classes in a bounding box.

    Features are queried from the public Overpass API. If `cache` is given,
    the result is stored there as JSON and re-used on later calls with the
    same region, tag and types, so the map does not change as OSM is edited
    and no network access is needed for re-plots.

    :type region: list of float
    :param region: [minlon, maxlon, minlat, maxlat] in decimal degrees
    :type tag: str
    :param tag: OSM tag key to filter on, e.g. "highway" or "waterway"
    :type types: list of str
    :param types: tag values to fetch, e.g. ["primary", "secondary"]
    :type cache: str
    :param cache: optional path to a JSON cache file
    :rtype: list of dict
    :return: one dictionary per OSM way with 'type' (the tag value), 'lon'
        and 'lat' (lists of vertex coordinates) keys
    """
    types = list(types)
    region = [float(v) for v in region[:4]]

    if cache and os.path.exists(cache):
        with open(cache) as f:
            cached = json.load(f)
        if (cached.get("region") == region and cached.get("tag") == tag
                and cached.get("types") == types):
            print(f"read {len(cached['lines'])} {tag} lines from {cache}")
            return cached["lines"]

    minlon, maxlon, minlat, maxlat = region
    query = (f'{OVERPASS_SETTINGS};'
             f'way["{tag}"~"^({"|".join(types)})$"]'
             f'({minlat},{minlon},{maxlat},{maxlon});'
             f'out geom;')
    print(f"fetching {tag} lines from OpenStreetMap ({', '.join(types)})")
    elements = overpass_query(query)

    lines = []
    for elem in elements:
        geometry = elem.get("geometry")
        if not geometry:
            continue
        lines.append({
            "type": elem.get("tags", {}).get(tag, "unknown"),
            "lon": [node["lon"] for node in geometry],
            "lat": [node["lat"] for node in geometry],
        })
    print(f"fetched {len(lines)} {tag} lines")

    if cache:
        with open(cache, "w") as f:
            json.dump({"region": region, "tag": tag, "types": types,
                       "lines": lines}, f)
    return lines


def fetch_osm_roads(
    region: List[float],
    road_types: Optional[List[str]] = None,
    cache: Optional[str] = None,
) -> List[Dict[str, object]]:
    """Fetch road lines inside a bounding box from OpenStreetMap.

    :type region: list of float
    :param region: [minlon, maxlon, minlat, maxlat] in decimal degrees
    :type road_types: list of str
    :param road_types: OSM 'highway' tag values to fetch, e.g. "primary",
        "residential", "track"; defaults to DEFAULT_ROAD_TYPES
    :type cache: str
    :param cache: optional path to a JSON cache file
    :rtype: list of dict
    :return: roads as returned by :func:`fetch_osm_lines`
    """
    return fetch_osm_lines(region, "highway",
                           road_types or DEFAULT_ROAD_TYPES, cache)


def fetch_osm_rivers(
    region: List[float],
    river_types: Optional[List[str]] = None,
    cache: Optional[str] = None,
) -> List[Dict[str, object]]:
    """Fetch river and stream lines inside a bounding box from OpenStreetMap.

    :type region: list of float
    :param region: [minlon, maxlon, minlat, maxlat] in decimal degrees
    :type river_types: list of str
    :param river_types: OSM 'waterway' tag values to fetch, e.g. "river",
        "stream", "canal"; defaults to DEFAULT_RIVER_TYPES
    :type cache: str
    :param cache: optional path to a JSON cache file
    :rtype: list of dict
    :return: rivers as returned by :func:`fetch_osm_lines`
    """
    return fetch_osm_lines(region, "waterway",
                           river_types or DEFAULT_RIVER_TYPES, cache)


def plot_line_groups(fig: pygmt.Figure, lines: List[Dict[str, object]],
                     rank: Dict[str, str], order: List[str],
                     pens: Dict[str, str]):
    """Draw OSM line features grouped by class, one plot call per group.

    Lines in the same group are joined into one multi-segment line (NaN rows
    separate segments). Groups are drawn in `order`, so later groups sit on
    top of earlier ones.

    :type fig: pygmt.Figure
    :param fig: figure with the map region and projection already set
    :type lines: list of dict
    :param lines: lines as returned by :func:`fetch_osm_lines`
    :type rank: dict
    :param rank: maps a line's 'type' to a group name; unlisted types fall
        in the first group of `order`
    :type order: list of str
    :param order: group names in drawing order
    :type pens: dict
    :param pens: pen for each group name
    """
    for group in order:
        xs, ys = [], []
        for line in lines:
            if rank.get(line["type"], order[0]) != group:
                continue
            xs += list(line["lon"]) + [float("nan")]
            ys += list(line["lat"]) + [float("nan")]
        if xs:
            fig.plot(x=xs, y=ys, pen=pens[group])


def plot_roads(fig: pygmt.Figure, roads: List[Dict[str, object]],
               road_styles: Dict[str, str]):
    """Draw road lines on a figure, with a pen chosen by road class.

    Minor roads are drawn first so larger roads sit on top.

    :type fig: pygmt.Figure
    :param fig: figure with the map region and projection already set
    :type roads: list of dict
    :param roads: roads as returned by :func:`fetch_osm_roads`
    :type road_styles: dict
    :param road_styles: 'major', 'secondary' and 'minor' pens; motorway,
        trunk and primary use 'major', secondary uses 'secondary', and all
        other classes use 'minor'
    """
    rank = {"motorway": "major", "motorway_link": "major",
            "trunk": "major", "trunk_link": "major",
            "primary": "major", "primary_link": "major",
            "secondary": "secondary", "secondary_link": "secondary"}
    plot_line_groups(fig, roads, rank, ["minor", "secondary", "major"],
                     road_styles)


def plot_rivers(fig: pygmt.Figure, rivers: List[Dict[str, object]],
                river_styles: Dict[str, str]):
    """Draw river lines on a figure, with a pen chosen by waterway class.

    :type fig: pygmt.Figure
    :param fig: figure with the map region and projection already set
    :type rivers: list of dict
    :param rivers: rivers as returned by :func:`fetch_osm_rivers`
    :type river_styles: dict
    :param river_styles: 'river' and 'minor' pens; rivers and canals use
        'river', and all other classes (streams, ditches...) use 'minor'
    """
    rank = {"river": "river", "canal": "river"}
    plot_line_groups(fig, rivers, rank, ["minor", "river"], river_styles)


def fetch_osm_coastline(
    region: List[float],
    cache: Optional[str] = None,
) -> List[Dict[str, object]]:
    """Fetch coastline lines inside a bounding box from OpenStreetMap.

    OSM coastline ways are drawn with land on the left and sea on the right
    of the direction of travel, which :func:`conform_grid_to_coastline`
    relies on.

    :type region: list of float
    :param region: [minlon, maxlon, minlat, maxlat] in decimal degrees
    :type cache: str
    :param cache: optional path to a JSON cache file
    :rtype: list of dict
    :return: coastline ways as returned by :func:`fetch_osm_lines`
    """
    return fetch_osm_lines(region, "natural", ["coastline"], cache)


def conform_grid_to_coastline(grid, coastline: List[Dict[str, object]],
                              band: float = 500.0):
    """Adjust a topography grid so its 0 m contour follows a coastline.

    GMT's high-resolution earth_relief grids blend land (SRTM) and ocean
    data, so their 0 m crossing can sit tens to hundreds of meters off the
    true shoreline. Within `band` meters of the coastline, cells on the land
    side with elevation <= 0 take the elevation of the nearest land cell (at
    least +1 m), and cells on the sea side with elevation >= 0 take the depth
    of the nearest sea cell (at most -1 m). Cells farther away are
    untouched, so the shading and contours elsewhere are unchanged.

    Land and sea are decided from the nearest coastline segment, using the
    OSM convention that land lies to the left of a coastline way.

    :type grid: xarray.DataArray
    :param grid: topography grid with (lat, lon) dimensions, as returned by
        pygmt.datasets.load_earth_relief
    :type coastline: list of dict
    :param coastline: coastline ways as returned by
        :func:`fetch_osm_coastline`
    :type band: float
    :param band: distance from the coastline in meters within which cells
        may be changed
    :rtype: xarray.DataArray
    :return: a copy of the grid with coastal cells adjusted
    """
    import numpy as np
    from scipy.ndimage import distance_transform_edt
    from scipy.spatial import cKDTree

    lat_name, lon_name = grid.dims
    lons = grid[lon_name].values
    lats = grid[lat_name].values

    # Local equirectangular projection to meters, fine at map scales
    mx = 111320.0 * np.cos(np.deg2rad(lats.mean()))
    my = 110540.0

    segs = []
    for way in coastline:
        x = np.asarray(way["lon"]) * mx
        y = np.asarray(way["lat"]) * my
        segs.append(np.column_stack([x[:-1], y[:-1], x[1:], y[1:]]))
    if not segs:
        return grid
    segs = np.vstack(segs)
    ax, ay, bx, by = segs.T
    dx, dy = bx - ax, by - ay
    seg_len2 = np.maximum(dx ** 2 + dy ** 2, 1e-12)
    half_len = 0.5 * np.sqrt(seg_len2.max())

    gx, gy = np.meshgrid(lons * mx, lats * my)
    px, py = gx.ravel(), gy.ravel()

    # Only cells whose nearest segment midpoint is close can be in the band
    tree = cKDTree(np.column_stack([(ax + bx) / 2, (ay + by) / 2]))
    d_mid, _ = tree.query(np.column_stack([px, py]), k=1)
    near = np.flatnonzero(d_mid < band + half_len)
    if near.size == 0:
        return grid

    # Exact point-to-segment distance over a few candidate segments
    k = min(8, len(segs))
    _, cand = tree.query(np.column_stack([px[near], py[near]]), k=k)
    cand = cand.reshape(len(near), k)
    qx, qy = px[near][:, None], py[near][:, None]
    t = np.clip(((qx - ax[cand]) * dx[cand] + (qy - ay[cand]) * dy[cand])
                / seg_len2[cand], 0, 1)
    dist = np.hypot(ax[cand] + t * dx[cand] - qx,
                    ay[cand] + t * dy[cand] - qy)
    dmin = dist.min(axis=1)

    # Side test: positive cross product means left of the segment (land).
    # When the closest point is a vertex shared by two segments, both are
    # equally close and may disagree, so sum their length-normalized cross
    # products (the angle-weighted normal at the vertex) to break the tie
    cross = (dx[cand] * (qy - ay[cand]) - dy[cand] * (qx - ax[cand]))
    ties = dist <= dmin[:, None] + 1e-6
    side = np.where(ties, cross / np.sqrt(seg_len2[cand]), 0).sum(axis=1)

    in_band = dmin < band
    land = side > 0

    z = grid.values.copy()
    flat = z.reshape(-1)
    to_land = np.zeros(flat.size, bool)
    to_sea = np.zeros(flat.size, bool)
    to_land[near[in_band & land]] = True
    to_sea[near[in_band & ~land]] = True
    # Exactly 0 m is ambiguous (SRTM marks water as 0), so flip it on either
    # side; otherwise the 0 m contour traces loops around flat 0 m patches
    to_land &= flat <= 0
    to_sea &= flat >= 0

    # Fill flipped cells from the nearest unflipped cell on the correct side
    # so colors and shading stay continuous, rather than a flat +/-1 m strip
    for flip, keep, clamp in ((to_land, flat > 0, 1.0),
                              (to_sea, flat < 0, -1.0)):
        if not flip.any() or not keep.any():
            flat[flip] = clamp
            continue
        source = (keep & ~to_land & ~to_sea).reshape(z.shape)
        _, (iy, ix) = distance_transform_edt(~source, return_indices=True)
        filled = z[iy, ix].reshape(-1)
        fill = filled[flip]
        flat[flip] = (np.maximum(fill, clamp) if clamp > 0
                      else np.minimum(fill, clamp))
    return grid.copy(data=z)


def get_default_styles() -> Dict[str, object]:
    """Get default style configuration for all map elements.

    Returns a dictionary containing all styling parameters for topography,
    coastlines, contours, markers, text, and bounding boxes. These can be
    overridden by passing a custom styles dict to plot_bbox_map().

    PyGMT Pen Specification Reference Guide.
    
    Overview:
        Pens control line appearance (thickness, color, style) for coastlines,
        contours, and other map features. Pen format: "width[unit],color[,style]"
    
    Components:
    -----------
    1. Width (line thickness):
        - Format: number + unit
        - Units: "p" (points, ~1/72 inch), "c" (cm), "i" (inches), "m" (mm)
        - Examples: "0.5p", "1p", "2p", "0.1c"
        - Common values:
          * Thin lines: "0.5p" or "1p"
          * Medium lines: "1.5p" or "2p"
          * Thick lines: "3p" or "4p"
    
    2. Color (line color):
        - Color names: "black", "red", "blue", "green", "white", etc.
        - Hex codes: "#FF0000" (red), "#00FF00" (green), "#0000FF" (blue)
        - RGB format: "255/0/0" for (R,G,B) values
        - Common colors: "gray", "lightgray", "darkgray", "yellow", "cyan", "magenta"
    
    3. Style (line pattern - OPTIONAL):
        - Solid (default): omit or use "-"
        - Dashed: use "-" with dash specification like "2p-" or "3p_5p" (3pt dash, 5pt gap)
        - Dotted: "1p.." (1pt dots)
        - Dash-dot: "2p-_1p" 
        - Note: Many GMT commands don't support style, only width and color
    
    PyGMT Pen Examples:
    -------------------
    Single parameter examples (most common):
        "0.5p,black"           → thin black line
        "1.5p,red"             → medium red line
        "2p,blue"              → thick blue line
        "1p,gray"              → thin gray line
        "0.75p,#FF6600"        → thin orange line (hex)
        "1.5p,100/150/255"     → medium custom-color line (RGB)
    
    With style (less common in PyGMT):
        "1p,black,-"           → 1pt black dashed line
        "1.5p,red,."           → 1.5pt red dotted line
    
    Typical Use in plot_bbox_map:
    ----------------------------
    bbox_pen="1.5p"            → 1.5 point thin line for bounding boxes
    bbox_pen="2p"              → 2 point medium line
    pen="0.5p,black"           → 0.5 point black line for contours
    shorelines=["1/0.8p,black"] → 1/0.8 point black shoreline
    
    Tips:
    -----
    - Start with "1p" or "1.5p" - adjust up or down for visibility
    - Use "black" or "gray" for professional maps
    - Avoid fancy styles; many PyGMT commands ignore them
    - Test on a small map area to see the effect before final plots

    Code	Style	        Format
    c	    Circle	        c0.5c (0.5 cm diameter)
    s	    Square	        s0.5c
    t	    Triangle	    t0.5c
    d	    Diamond	        d0.5c
    h	    Hexagon	        h0.5c
    p	    Pentagon	    p0.5c
    x	    X/cross	        x0.5c
    +	    Plus sign	    +0.5c
    *	    Asterisk	    *0.5c
    a	    Star (5-point)	a0.5c

    Returns
    -------
    dict
        Dictionary with keys for each map element's styling:
        - topo: topography image settings
        - contour: contour line settings
        - coastline: coastline settings
        - marker: default marker settings
        - text: default text label settings
        - bbox: bounding box settings
    """
    return {
        # Topography (grdimage) settings
        "topo": {
            "cmap": "gmt/geo",
            "shading": True,
        },
        # Contour (grdcontour) settings
        "contour": {
            "pen": "0.5p,black",
            "annotate_size": "8p",  # Font size for contour labels
        },
        # Coastline (coast) settings
        "coastline": {
            "shoreline_width": "0.8p",
            "shoreline_color": "black",
        },
        # Default marker (plot) settings
        "marker": {
            "style": "c0.4c",  # Circle, 0.4 cm diameter
            "color": "red",
            "pen": "0.5p,black",
        },
        # Default text label settings
        "text": {
            "font": "8p,Helvetica",  # 8 point Helvetica (color set separately)
            "color": "black",        # Text label color
            "offset": "0.15c/0.15c",  # Offset from marker point
        },
        # Bounding box settings
        "bbox": {
            "color": "red",
            "pen_width": "1.5p",
        },
        # Road pens by class (see plot_roads)
        "road": {
            "major": "1.25p,30/30/30",
            "secondary": "0.9p,30/30/30",
            "minor": "0.6p,60/60/60",
        },
        # River pens by class (see plot_rivers)
        "river": {
            "river": "0.2p,30/100/200",
            "minor": "0.1p,30/100/200",
        },
    }


def plot_bbox_map(
    region: List[float] | List[List[float]],
    output: str = "bbox_map.png",
    topo_resolution: str = "01m",
    projection: str = "M12c",
    contour_interval: float = 500,
    map_scale: str = None,
    markers: Optional[List[Dict[str, object]]] = None,
    title: Optional[str] = None,
    show: bool = True,
    utm_zone: Optional[int] = None,
    utm_hemisphere: str = "N",
    topo_data_source: str = "gebco",
    coast_resolution: str = "i",
    river_resolution: Optional[str] = None,
    border_resolution: Optional[str] = None,
    styles: Optional[Dict[str, object]] = None,
    coastline_from_grid: bool = False,
    colorbar: bool = False,
    geo_format: Optional[str] = None,
    roads: Optional[List[Dict[str, object]]] = None,
    rivers: Optional[List[Dict[str, object]]] = None,
    coastline: Optional[List[Dict[str, object]]] = None,
    coast_band: float = 500.0,
) -> pygmt.Figure:
    """Plot a bounding-box map with topography, contours, coastlines, and markers.

    Parameters
    ----------
    region
        Either a single region [minlon, maxlon, minlat, maxlat] or a list of regions.
        If multiple regions are provided, the largest becomes the main map and
        smaller ones are drawn as outlined boxes on top.
    output
        Path to save the figure.
    topo_resolution
        GMT earth_relief resolution string, e.g. "01m", "03m", "06m".
    projection
        PyGMT projection string, such as "M12c" for Mercator 12 cm width.
    cmap
        Colormap for topography.
    contour_interval
        Elevation contour interval in meters.
    map_scale
        Optional map scale string for basemap, e.g. "jBL+w100k+f+o0.5c/0.5c".
    markers
        Optional list of marker dictionaries, each with keys:
        - lon/lat: longitude/latitude (lon/lat coordinates)
        - easting/northing: UTM easting/northing (if utm_zone is specified)
        - style: GMT symbol style string (default from styles["marker"]["style"])
        - color: fill color (default from styles["marker"]["color"])
        - pen: outline pen (default from styles["marker"]["pen"])
        - label: optional text label
        - text_color: text label color (default from styles["text"]["color"], overrides per-marker)
        - label_offset: text label offset (default from styles["text"]["offset"])
        - font: text font specification (default from styles["text"]["font"])
    title
        Optional title text.
    show
        If True, show the figure after saving.
    utm_zone
        UTM zone number (1-60). If provided, region and markers are assumed
        to be in UTM coordinates and will be converted to lon/lat.
    utm_hemisphere
        "N" for Northern hemisphere, "S" for Southern hemisphere.
        Only used if utm_zone is specified.
    topo_data_source
        Data source for topography, e.g. "gebco" (default) or "etopo1g".
    coast_resolution
        GMT coastline resolution. Options from highest to lowest detail:
        "f" (full), "h" (high), "i" (intermediate, default), "l" (low), "c" (crude).
    river_resolution
        GMT river resolution. Same resolution levels as coast_resolution.
        If None, rivers are not displayed.
    border_resolution
        GMT political boundary resolution. Same resolution levels as coast_resolution.
        If None, borders are not displayed.
    bbox_color
        Color for the outlined bounding boxes (default "red").
        Note: Overridden by styles["bbox"]["color"] if provided.
    bbox_pen
        Pen specification for bounding boxes (default "1.5p" = 1.5 point line).
        Note: Overridden by styles["bbox"]["pen_width"] if provided.
    styles
        Optional dictionary of style overrides. Keys can include "topo", "contour",
        "coastline", "marker", "text", "bbox". Unspecified values use defaults from
        get_default_styles(). Example:
            styles = {
                "contour": {"pen": "1p,blue", "annotate_size": "10p"},
                "marker": {"style": "s0.5c", "color": "green"},
            }
    coastline_from_grid
        If True, draw the coastline as the 0 m contour of the topography grid
        instead of the GSHHG shoreline database, so the coastline matches the
        land/sea boundary of the plotted topography exactly.
    colorbar
        If True, draw a vertical elevation colorbar to the right of the map.
    geo_format
        Optional GMT FORMAT_GEO_MAP string for the axis annotations, e.g.
        "ddd.xxF" for decimal degrees with a hemisphere letter (155.30°W).
        If None, GMT's default degree-minute format is used.
    roads
        Optional road lines from :func:`fetch_osm_roads`, drawn above the
        topography and below the markers. Pens come from styles["road"].
    rivers
        Optional river lines from :func:`fetch_osm_rivers`, drawn below the
        roads. Pens come from styles["river"]. For GMT's built-in (coarser)
        GSHHG rivers use `river_resolution` instead.
    coastline
        Optional coastline ways from :func:`fetch_osm_coastline`. If given,
        the topography is adjusted near the coast with
        :func:`conform_grid_to_coastline` so the land/sea colors and the
        coastline (drawn as the 0 m contour) match the OSM shoreline, and
        therefore OSM roads and rivers. Implies `coastline_from_grid`.
    coast_band
        Distance in meters from the coastline within which the topography
        may be adjusted when `coastline` is given (default 500).

    Returns
    -------
    pygmt.Figure
        The generated figure.
    """
    # Get default styles and merge with user-provided overrides
    default_styles = get_default_styles()
    if styles is None:
        styles = {}
    
    
    # Deep merge: update defaults with any provided overrides
    merged_styles = {}
    for key in default_styles:
        merged_styles[key] = {**default_styles[key], **styles.get(key, {})}
    styles = merged_styles

    # Handle multiple regions: main region should be the first
    if region and isinstance(region[0], (list, tuple)):
        # Multiple regions provided
        main_region = region[0]
        smaller_regions = region[1:]
    else:
        # Single region provided
        main_region = region
        smaller_regions = []

    # Convert region from UTM to lon/lat if needed
    if detect_coordinate_system(main_region) == "utm":
        region_lonlat = convert_region_to_lonlat(main_region, utm_zone, utm_hemisphere)
    else:
        region_lonlat = main_region

    # Convert smaller regions from UTM to lon/lat if needed
    smaller_regions_lonlat = []
    if smaller_regions:
        for small_region in smaller_regions:
            if detect_coordinate_system(small_region) == "utm":
                small_region_lonlat = convert_region_to_lonlat(small_region, utm_zone, utm_hemisphere)
                smaller_regions_lonlat.append(small_region_lonlat)
            else:
                smaller_regions_lonlat.append(small_region)

    # Convert markers from UTM to lon/lat if needed
    markers_lonlat = markers
    if markers is not None and utm_zone is not None:
        markers_lonlat = convert_markers_to_lonlat(markers, utm_zone, utm_hemisphere)

    fig = pygmt.Figure()

    # Load satellite-derived topography/bathymetry from the GMT remote repository.
    topo_grid = pygmt.datasets.load_earth_relief(
        resolution=topo_resolution,
        region=region_lonlat,
        data_source=topo_data_source,
    )

    # Snap the DEM's land/sea boundary to the vector coastline
    if coastline:
        topo_grid = conform_grid_to_coastline(topo_grid, coastline,
                                              coast_band)
        coastline_from_grid = True

    # Map frame is drawn here only, so the annotation format applies here.
    # GMT's automatic ticks fall on arc-minutes, which give awkward decimal
    # labels (e.g. 155.33), so pick round decimal-degree spacings instead
    if geo_format:
        geo_config = pygmt.config(FORMAT_GEO_MAP=geo_format)
        dx = nice_interval(region_lonlat[1] - region_lonlat[0])
        dy = nice_interval(region_lonlat[3] - region_lonlat[2])
        frame = [f"xa{dx}f{dx / 2}", f"ya{dy}f{dy / 2}"]
    else:
        geo_config = contextlib.nullcontext()
        frame = ["a"]
    if title:
        frame.append(f"+t{title}")
    with geo_config:
        fig.grdimage(
            grid=topo_grid,
            projection=projection,
            region=region_lonlat,
            cmap="gmt/geo",
            shading=styles["topo"]["shading"],
            frame=frame,
        )

    if colorbar:
        fig.colorbar(position="JMR+o0.5c/0c+v",
                     frame=["xaf+lElevation", "y+lm"])

    # Contour the topography (excluding sea level at 0 meters to avoid overlap with coastline)
    # Generate contour levels, excluding 0
    contour_levels = []
    # Negative levels (below sea level)
    level = -6000
    while level < -contour_interval/2:
        contour_levels.append(level)
        level += contour_interval
    # Positive levels (above sea level)
    level = contour_interval
    while level <= 9000:
        contour_levels.append(level)
        level += contour_interval
    
    # Select every other contour level for annotation to avoid overcrowding
    annotate_levels = contour_levels[::2]
    
    fig.grdcontour(
        grid=topo_grid,
        levels=contour_levels,
        pen=styles["contour"]["pen"],
        annotation=annotate_levels,
    )

    coast_pen = (f"{styles['coastline']['shoreline_width']},"
                 f"{styles['coastline']['shoreline_color']}")

    # Coastline as the 0 m contour of the same grid used for shading
    if coastline_from_grid:
        fig.grdcontour(grid=topo_grid, levels="0,", pen=coast_pen)

    # Add coastlines over the topo layer.
    coast_args = {
        "region": region_lonlat,
        "projection": projection,
        "resolution": coast_resolution,
    }
    if not coastline_from_grid:
        coast_args["shorelines"] = [f"1/{coast_pen}"]

    # Add rivers if specified
    if river_resolution:
        coast_args["rivers"] = river_resolution

    # Add political boundaries if specified
    if border_resolution:
        coast_args["borders"] = border_resolution

    # coast() needs at least one feature to draw
    if any(k in coast_args for k in ("shorelines", "rivers", "borders")):
        fig.coast(**coast_args)

    # Rivers below roads so road crossings sit on top
    if rivers:
        plot_rivers(fig, rivers, styles["river"])
    if roads:
        plot_roads(fig, roads, styles["road"])

    # Scale Bar
    if map_scale:
        fig.basemap(map_scale=map_scale)

    # Draw outlined boxes for smaller regions
    if smaller_regions_lonlat:
        for small_region in smaller_regions_lonlat:
            minlon, maxlon, minlat, maxlat = small_region[:4]
            # Create rectangle coordinates: west, east, south, north
            fig.plot(
                region=region_lonlat,
                projection=projection,
                x=[minlon, maxlon, maxlon, minlon, minlon],
                y=[minlat, minlat, maxlat, maxlat, minlat],
                pen=f"{styles['bbox']['pen_width']},{styles['bbox']['color']}",
            )

    # Plot markers and optional labels
    if markers_lonlat:
        for marker in markers_lonlat:
            lon = float(marker["lon"])
            lat = float(marker["lat"])
            style = marker.get("style", styles["marker"]["style"])
            fill = marker.get("color", styles["marker"]["color"])
            pen = marker.get("pen", styles["marker"]["pen"])
            fig.plot(x=lon, y=lat, style=style, fill=fill, pen=pen)

            label = marker.get("label")
            if label:
                offset = marker.get("label_offset", styles["text"]["offset"])
                font = marker.get("font", styles["text"]["font"])
                fig.text(x=lon, y=lat, text=label, font=font, offset=offset)

    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    fig.savefig(output)
    if show:
        fig.show(method="external")
    return fig
    

def simblast_regional(main_region):
    """Example usage of the plot_bbox_map function with multiple bounding boxes."""
    # Define main region (largest) and two smaller regions
    # main_region = [115.0, 140.5, 33.5, 45.5]
    high_res_region = [425_000.000, 565_000.000, 4_490_000.000, 4_680_000.000]

    # Convert all to Lon/Lat if needed
    regions = [main_region, high_res_region]


    station_color = "red"
    city_color = "orange"
    markers_lonlat_example = [
        {"label": "NKNTS", "lon": 129.0297, "lat": 41.33, 
         "font": "10p,Helvetica,white", "color": "yellow", 
         "pen": "1p,black", "style": "a0.5c"},         
         # STATIONS
        {"label": "IC.MDJ", "lon": 129.5934, "lat": 44.6176, 
          "font": "10p,Helvetica,white", "color": station_color, 
          "pen": "1p,black", "style": "t0.5c"},
        {"label": "IU.INCN", "lon": 126.6244, "lat": 37.4777, 
          "font": "10p,Helvetica,white", "color": station_color, 
          "pen": "1p,black", "style": "t0.5c"},
        {"label": "IC.BJT", "lon": 116.1679, "lat": 40.0183, 
          "font": "10p,Helvetica,white", "color": station_color, 
          "pen": "1p,black", "style": "t0.5c"},
        {"label": "G.INU", "lon": 137.029, "lat": 35.35, 
          "font": "10p,Helvetica,white", "color": station_color, 
          "pen": "1p,black", "style": "t0.5c"},
        # CITIES
        # {"label": "Seoul", "lon": 126.978, "lat": 37.5665, 
        #   "font": "10p,Helvetica,white", "color": station_color, 
        #   "pen": "1p,black", "style": "c0.5c"},
        # {"label": "Tokyo", "lon": 139.6917, "lat": 35.6895, 
        #   "font": "10p,Helvetica,white", "color": city_color, 
        #   "pen": "1p,black", "style": "c0.5c"}
    ]

    plot_bbox_map(
        region=regions,
        output="simblast_regional.png",
        topo_resolution="01m",
        projection="M12c",
        contour_interval=10000,
        map_scale="jBR+k100k+f+o0.5c/0.5c",
        utm_zone=52,
        utm_hemisphere="N",
        markers=markers_lonlat_example,
        title=None,
        coast_resolution="h",
        styles={},  # Pass all styles at once
    )


def simblast_high_res():
    """Example usage of the plot_bbox_map function with multiple bounding boxes."""
    region = [425_000.000, 565_000.000, 4_490_000.000, 4_680_000.000]

    markers_lonlat_example = [
        {"label": "NKNTS", "lon": 129.0297, "lat": 41.33, 
         "font": "10p,Helvetica,red", "color": "yellow", 
         "pen": "1p,black", "style": "a0.5c"},        
    ]

    plot_bbox_map(
        region=region,
        output="simblast_high_res.png",
        topo_resolution="01m",
        projection="M12c",
        contour_interval=250,
        map_scale="jBL+w100k+f+o0.5c/0.5c",
        utm_zone=52,
        utm_hemisphere="N",
        markers=markers_lonlat_example,
        title=None,
        coast_resolution="h",
        styles={},  # Pass all styles at once
    )


def gulkana():
    """Example usage of the plot_bbox_map function with multiple bounding boxes."""
    region = [-145.59971,-145.3003,63.198378,63.349744]

    markers_lonlat_example = [
        {"label": "NKNTS", "lon": 129.0297, "lat": 41.33, 
         "font": "10p,Helvetica,red", "color": "yellow", 
         "pen": "1p,black", "style": "a0.5c"},        
    ]

    plot_bbox_map(
        region=region,
        output="cryo.png",
        topo_data_source="igpp",
        topo_resolution="01m",
        projection="M12c",
        contour_interval=250,
        map_scale="jBL+w10k+f+o0.5c/0.5c",
        utm_zone=6,
        utm_hemisphere="N",
        title=None,
        coast_resolution="h",
        styles={},  # Pass all styles at once
    )

def read_kmz(path: str) -> List[Dict[str, object]]:
    """Read point Placemarks from a KMZ (or plain KML) file.

    Each Placemark that contains a Point is returned as a marker dictionary
    with 'lon', 'lat' and 'label' (the Placemark name) keys, ready to be
    passed to :func:`plot_bbox_map`. Non-point geometries are skipped.

    .. note::
        Not used by :func:`main`, which reads JSON points files; kept for
        posterity, e.g. for converting KMZ station files to JSON.

    :type path: str
    :param path: path to a .kmz archive or a .kml file
    :rtype: list of dict
    :return: marker dictionaries, one per point Placemark
    """
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as kmz:
            kml_names = [n for n in kmz.namelist()
                         if n.lower().endswith(".kml")]
            if not kml_names:
                raise ValueError(f"no .kml file found inside {path}")
            # Google Earth convention is a root-level doc.kml
            kml_name = "doc.kml" if "doc.kml" in kml_names else kml_names[0]
            root = ET.fromstring(kmz.read(kml_name))
    else:
        root = ET.parse(path).getroot()

    def _local(tag):
        """Strip the XML namespace from a tag name."""
        return tag.rsplit("}", 1)[-1]

    markers = []
    for elem in root.iter():
        if _local(elem.tag) != "Placemark":
            continue
        name, coords = None, None
        for child in elem.iter():
            tag = _local(child.tag)
            if tag == "name" and name is None:
                name = (child.text or "").strip()
            elif tag == "Point":
                for sub in child.iter():
                    if _local(sub.tag) == "coordinates" and sub.text:
                        coords = sub.text.strip().split()[0]
        if coords is None:
            continue
        lon, lat = [float(v) for v in coords.split(",")[:2]]
        markers.append({"lon": lon, "lat": lat, "label": name})

    if not markers:
        raise ValueError(f"no point Placemarks found in {path}")
    return markers


# Readable marker names mapped to GMT symbol codes
MARKER_CODES = {
    "circle": "c", "square": "s", "triangle": "t",
    "inverted_triangle": "i", "diamond": "d", "star": "a",
    "hexagon": "h", "pentagon": "n", "octagon": "g",
    "cross": "x", "plus": "+",
}

# Points always plotted by main(), before any points from JSON files. Each
# entry is a dictionary; only 'lat' and 'lon' are required. See
# make_point_marker() for all keys and defaults. Example:
#     {"lat": 19.4069, "lon": -155.2834, "label": "HVO",
#      "color": "yellow", "marker": "star", "size": "0.5c",
#      "text_color": "white", "font_size": "10p"}
POINTS = [
]


def make_point_marker(point: Dict[str, object]) -> Dict[str, object]:
    """Convert a user point dictionary into a plot_bbox_map marker.

    Recognized keys (only 'lat' and 'lon' are required):

    - lat, lon: coordinates in decimal degrees
    - label: text label drawn next to the point (default: no label)
    - color: marker fill color (default "yellow")
    - marker: a name from MARKER_CODES (e.g. "star", "circle") or a raw GMT
      symbol code (e.g. "a") (default "circle")
    - size: marker size with unit (default "0.4c")
    - pen: marker outline pen (default "0.75p,black")
    - text_color: label color (default "black")
    - font_size: label font size (default "8p")
    - font: label font name (default "Helvetica")
    - label_offset: label offset from the point (default "0.3c/0.3c")

    :type point: dict
    :param point: user point dictionary as described above
    :rtype: dict
    :return: marker dictionary accepted by :func:`plot_bbox_map`
    """
    if "lat" not in point or "lon" not in point:
        raise ValueError(f"point is missing 'lat' or 'lon': {point}")
    marker = point.get("marker", "circle")
    code = MARKER_CODES.get(marker, marker)
    font = (f"{point.get('font_size', '8p')},"
            f"{point.get('font', 'Helvetica')},"
            f"{point.get('text_color', 'black')}")
    return {
        "lon": float(point["lon"]),
        "lat": float(point["lat"]),
        "label": point.get("label"),
        "style": f"{code}{point.get('size', '0.4c')}",
        "color": point.get("color", "yellow"),
        "pen": point.get("pen", "0.75p,black"),
        "font": font,
        "label_offset": point.get("label_offset", "0.3c/0.3c"),
    }


def recolor_pens(pens: Dict[str, str], color: str) -> Dict[str, str]:
    """Replace the color of every pen in a dictionary, keeping widths.

    :type pens: dict
    :param pens: pens like {"major": "1.25p,30/30/30", ...}
    :type color: str
    :param color: any GMT color, e.g. a name ("white", "steelblue"), an
        R/G/B triplet ("200/0/0") or a hex code ("#FF8800")
    :rtype: dict
    :return: the same keys with pens of the form "<width>,<color>"
    """
    return {key: f"{pen.split(',')[0]},{color}" for key, pen in pens.items()}


def main(points: Optional[List[Dict[str, object]]] = None):
    """Plot a map of points from one or more JSON points files.

    Each JSON file holds a list of point dictionaries (see
    :func:`make_point_marker`). Points are drawn on topography with
    :func:`plot_bbox_map` in the order the files are given, so later files
    plot on top. The map region must be given explicitly with --region.

    :type points: list of dict
    :param points: points plotted before those from the JSON files, each a
        dictionary as described in :func:`make_point_marker`; defaults to the
        module-level POINTS list
    """
    parser = argparse.ArgumentParser(
        description="Plot a map of points from JSON points files")
    parser.add_argument("points_files", nargs="+", metavar="POINTS_JSON",
                        help="JSON file(s) with a list of point dictionaries")
    parser.add_argument("-r", "--region", nargs=4, type=float, required=True,
                        metavar=("MINLON", "MAXLON", "MINLAT", "MAXLAT"),
                        help="map extent in decimal degrees")
    parser.add_argument("-o", "--output", default=None,
                        help="output figure path (default: <first points "
                             "file name>.png)")
    parser.add_argument("--topo_resolution", default="15s",
                        help="earth_relief resolution, e.g. 01s, 03s, 15s "
                             "(default 15s)")
    parser.add_argument("--topo_data_source", default=None,
                        choices=["igpp", "gebco", "gebcosi", "synbath"],
                        help="relief data source; default is igpp (SRTM on "
                             "land) for 01s/03s and gebco otherwise")
    parser.add_argument("--coast_source", default="grid",
                        choices=["grid", "gshhg", "osm"],
                        help="draw coastline as the 0 m contour of the "
                             "topography (grid, default), from the GSHHG "
                             "shoreline database (gshhg), or from the "
                             "OpenStreetMap coastline with the topography "
                             "adjusted near the coast to match it (osm; use "
                             "this with --roads/--rivers so they line up)")
    parser.add_argument("--coast_cache", default=None,
                        help="JSON file caching the OSM coastline (default: "
                             "<output name>_coast.json); delete it to "
                             "re-download")
    parser.add_argument("--coast_band", type=float, default=500.0,
                        help="with --coast_source osm, distance in meters "
                             "from the coast within which the topography "
                             "may be adjusted (default 500)")
    parser.add_argument("--contour_interval", type=float, default=500,
                        help="contour interval in meters (default 500)")
    parser.add_argument("--map_scale", default="jBL+w10k+f+u+o0.5c/0.5c",
                        help="GMT map scale string; +u appends units to "
                             "the labels")
    parser.add_argument("--no_colorbar", action="store_true",
                        help="do not draw the elevation colorbar")
    parser.add_argument("--geo_format", default="ddd.xxF",
                        help="GMT FORMAT_GEO_MAP for axis labels (default "
                             "ddd.xxF, decimal degrees; ddd:mmF for "
                             "degree-minutes)")
    parser.add_argument("--projection", default="M15c",
                        help="GMT projection (default M15c)")
    parser.add_argument("--title", default=None, help="map title")
    parser.add_argument("--no_labels", action="store_true",
                        help="do not draw point labels")
    parser.add_argument("--roads", action="store_true",
                        help="draw OpenStreetMap roads inside the region")
    parser.add_argument("--road_types", nargs="+", default=None,
                        metavar="TYPE",
                        help="OSM highway classes to draw (default: "
                             f"{' '.join(DEFAULT_ROAD_TYPES)}); add e.g. "
                             "residential, unclassified or track for "
                             "smaller roads")
    parser.add_argument("--roads_cache", default=None,
                        help="JSON file caching the fetched roads (default: "
                             "<output name>_roads.json); delete it to "
                             "re-download")
    parser.add_argument("--road_color", default=None,
                        help="road color, e.g. white, black, orange, "
                             "darkred, or R/G/B or #RRGGBB (default dark "
                             "gray); line widths are unchanged")
    parser.add_argument("--rivers", action="store_true",
                        help="draw OpenStreetMap rivers inside the region")
    parser.add_argument("--river_types", nargs="+", default=None,
                        metavar="TYPE",
                        help="OSM waterway classes to draw (default: "
                             f"{' '.join(DEFAULT_RIVER_TYPES)}); add e.g. "
                             "stream for smaller channels")
    parser.add_argument("--rivers_cache", default=None,
                        help="JSON file caching the fetched rivers (default: "
                             "<output name>_rivers.json); delete it to "
                             "re-download")
    parser.add_argument("--river_color", default=None,
                        help="river color, e.g. blue, dodgerblue, cyan, "
                             "navy, or R/G/B or #RRGGBB (default medium "
                             "blue); line widths are unchanged")
    parser.add_argument("--show", action="store_true",
                        help="open the figure after saving")
    args = parser.parse_args()

    points = list(POINTS if points is None else points)
    for fid in args.points_files:
        with open(fid) as f:
            file_points = json.load(f)
        print(f"read {len(file_points)} points from {fid}")
        points += file_points
    markers = [make_point_marker(p) for p in points]
    if args.no_labels:
        for marker in markers:
            marker["label"] = None

    region = args.region
    output = (args.output
              or os.path.splitext(args.points_files[0])[0] + ".png")

    # GEBCO stops at 15s; only igpp (SRTM-backed) provides 01s and 03s
    topo_data_source = args.topo_data_source
    if topo_data_source is None:
        topo_data_source = ("igpp" if args.topo_resolution in ("01s", "03s")
                            else "gebco")
    print(f"topography: {topo_data_source} @ {args.topo_resolution}")

    roads = None
    if args.roads:
        cache = (args.roads_cache
                 or os.path.splitext(output)[0] + "_roads.json")
        roads = fetch_osm_roads(region, args.road_types, cache=cache)

    rivers = None
    if args.rivers:
        cache = (args.rivers_cache
                 or os.path.splitext(output)[0] + "_rivers.json")
        rivers = fetch_osm_rivers(region, args.river_types, cache=cache)

    coastline = None
    if args.coast_source == "osm":
        cache = (args.coast_cache
                 or os.path.splitext(output)[0] + "_coast.json")
        coastline = fetch_osm_coastline(region, cache=cache)
        if not coastline:
            print("no OSM coastline in region; using the topography's "
                  "0 m contour")

    styles = {}
    defaults = get_default_styles()
    if args.road_color:
        styles["road"] = recolor_pens(defaults["road"], args.road_color)
    if args.river_color:
        styles["river"] = recolor_pens(defaults["river"], args.river_color)

    plot_bbox_map(
        region=region,
        output=output,
        topo_resolution=args.topo_resolution,
        topo_data_source=topo_data_source,
        projection=args.projection,
        contour_interval=args.contour_interval,
        map_scale=args.map_scale,
        markers=markers,
        title=args.title,
        show=args.show,
        coast_resolution="f",
        coastline_from_grid=args.coast_source in ("grid", "osm"),
        colorbar=not args.no_colorbar,
        geo_format=args.geo_format,
        roads=roads,
        rivers=rivers,
        styles=styles,
        coastline=coastline,
        coast_band=args.coast_band,
    )
    print(f"saved map to {output}")


if __name__ == "__main__":
    main()
