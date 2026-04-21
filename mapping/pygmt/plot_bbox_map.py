"""Create a PyGMT map for a bounding box and add markers for stations/cities.

This script plots topography, contour elevations, coastlines, and allows
simple latitude/longitude markers with optional labels. Supports both lon/lat
and UTM coordinate systems.

* Vibecoded by GitHub Copilot
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

    fig.grdimage(
        grid=topo_grid,
        projection=projection,
        region=region_lonlat,
        cmap="gmt/geo",
        shading=styles["topo"]["shading"],
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
        pen=styles["contour"]["pen"],
        annotation=annotate_levels,
    )

    # Add coastlines over the topo layer.
    coast_args = {
        "region": region_lonlat,
        "projection": projection,
        "shorelines": [f"1/{styles['coastline']['shoreline_width']},{styles['coastline']['shoreline_color']}"],
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
        map_scale="jBR+w100k+f+o0.5c/0.5c",
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


if __name__ == "__main__":
    INCN = [126.0, 130., 37.25, 42.5]
    MDJ = [128.0, 131., 40.5, 45.0]
    simblast_regional(MDJ)
    # simblast_high_res()
