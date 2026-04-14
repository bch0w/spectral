"""Create a PyGMT map for a bounding box and add markers for stations/cities.

This script plots topography, contour elevations, coastlines, and allows
simple latitude/longitude markers with optional labels. Supports both lon/lat
and UTM coordinate systems.

Vibecoded by GitHub Copilot
"""
import os
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


def plot_bbox_map(
    region: List[float],
    output: str = "bbox_map.png",
    topo_resolution: str = "01m",
    projection: str = "M12c",
    cmap: str = "geo",
    contour_interval: float = 500,
    markers: Optional[List[Dict[str, object]]] = None,
    title: Optional[str] = None,
    show: bool = True,
    utm_zone: Optional[int] = None,
    utm_hemisphere: str = "N",
    topo_data_source: str = "gebco",
    coast_resolution: str = "i",
    river_resolution: Optional[str] = None,
    border_resolution: Optional[str] = None,
) -> pygmt.Figure:
    """Plot a bounding-box map with topography, contours, coastlines, and markers.

    Parameters
    ----------
    region
        [minlon, maxlon, minlat, maxlat] in lon/lat coordinates, or
        [min_easting, max_easting, min_northing, max_northing] if utm_zone is specified.
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
    markers
        Optional list of marker dictionaries, each with keys:
        - lon/lat: longitude/latitude (lon/lat coordinates)
        - easting/northing: UTM easting/northing (if utm_zone is specified)
        - style: GMT symbol style string (default "c0.4c")
        - color: fill color
        - pen: outline pen
        - label: optional text label
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

    Returns
    -------
    pygmt.Figure
        The generated figure.
    """
    # Convert region from UTM to lon/lat if needed
    region_lonlat = convert_region_to_lonlat(region, utm_zone, utm_hemisphere)

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

    fig.grdimage(
        grid=topo_grid,
        projection=projection,
        region=region_lonlat,
        cmap="gmt/geo",
        shading=True,
        frame=["a"] + ([f"+t{title}"] if title else []),
    )

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
        pen="0.5p,black",
        annotation=annotate_levels,
    )

    # Add coastlines over the topo layer.
    coast_args = {
        "region": region_lonlat,
        "projection": projection,
        "shorelines": ["1/0.8p,black"],
        "frame": ["a"] + ([f"+t{title}"] if title else []),
        "resolution": coast_resolution,
    }
    
    # Add rivers if specified
    if river_resolution:
        coast_args["rivers"] = river_resolution

    # Add political boundaries if specified
    if border_resolution:
        coast_args["borders"] = border_resolution

    fig.coast(**coast_args)

    # Plot markers and optional labels
    if markers_lonlat:
        for marker in markers_lonlat:
            lon = float(marker["lon"])
            lat = float(marker["lat"])
            style = marker.get("style", "c0.4c")
            fill = marker.get("color", "red")
            pen = marker.get("pen", "0.5p,black")
            fig.plot(x=lon, y=lat, style=style, fill=fill, pen=pen)

            label = marker.get("label")
            if label:
                offset = marker.get("label_offset", "0.15c")
                font = marker.get("font", "8p,black,Helvetica-Bold")
                fig.text(x=lon, y=lat, text=label, font=font, offset=f"0.15c/0.15c")

    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    fig.savefig(output)
    if show:
        fig.show(method="external")
    return fig
    
def example_usage():
    """Example usage of the plot_bbox_map function."""
    region_lonlat_example = [125.5, 130.5, 34, 42.5]
    markers_lonlat_example = [
        {"lon": 129.0297, "lat": 41.33, "label": "NKNTS", "color": "red"},
    ]

    plot_bbox_map(
        region=region_lonlat_example,
        output="mapping/pygmt/simblast.png",
        topo_resolution="01m",
        projection="M12c",
        contour_interval=200,
        markers=markers_lonlat_example,
        title="SimBlast NK",
        coast_resolution="h",  # Use high resolution coastlines
    )

def simblast_paper():
    """Example with UTM coordinates (Zone 11N)
    Same location as above, converted to UTM Zone 11N
    """
    region_utm_example = [425_000.000, 565_000.000, 4_490_000.000, 4_680_000.000]
    markers_utm_example = [
        {"lon": 129.0297, "lat": 41.33, "label": "NKNTS", "color": "red"},
    ]

    plot_bbox_map(
        region=region_utm_example,
        output="simblast_utm.png",
        topo_resolution="01m",
        projection="M12c",
        contour_interval=250,
        markers=markers_utm_example,
        title=None,
        utm_zone=52,
        utm_hemisphere="N",
        coast_resolution="h",  # Use high resolution coastlines
    )


if __name__ == "__main__":
    simblast_paper()
