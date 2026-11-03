"""
Inset map showing a globe with given lat/lon locations to highlight a region
of interest from a global perspective.

Points, lines of constant latitude, and map styling are read from a JSON
file. Every key except "points" is optional and falls back to `DEFAULTS`.
Example (see globe_inset_example.json):

    {
        "center": [-145, 63],
        "width": "12c",
        "land": "seashell",
        "water": "lightblue",
        "parallels": [{"lat": 66.5607, "pen": "1p,purple,--"}],
        "points": [
            {"lat": 63.27, "lon": -145.45, "style": "d0.2c",
             "fill": "red", "pen": "black", "label": "Gulkana Glacier"}
        ]
    }

Point keys: lat, lon (required); style, fill, pen, label, font, justify,
offset, arrow, arrow_pen, arrow_fill, arrow_head, arrow_gap (optional, see
`DEFAULT_POINT`). Set "arrow": true to draw an arrow from an offset label
to its marker. If "center" is omitted the map is centered on the mean
location of the points and catalog events.

"projection" picks the map projection by name: orthographic (default),
azimuthal_equidistant, lambert_azimuthal, robinson, mollweide, hammer,
winkel_tripel, eckert4, eckert6, equidistant_cylindrical, miller or
mercator. Any other value is used as a full GMT -J string, e.g.
"G-150/45/10c", in which case "width" is ignored. "region" overrides
the default whole-world GMT -R region.

Earthquake catalogs are read from QuakeML with ObsPy and listed under
"catalogs". Paths may be globs, and relative paths are relative to the
JSON file:

    "catalogs": [
        {"files": ["catalogs/*.xml"], "symbol": "c",
         "color_by": "depth", "cmap": "viridis", "cmap_reverse": true,
         "size_by_magnitude": true, "size_range": [0.05, 0.3]}
    ]

Catalog keys: files (required); symbol, size, size_by_magnitude,
size_range, mag_range, color_by ("depth", "time" or null), fill, cmap,
cmap_reverse, color_range, colorbar, colorbar_label, pen, name, legend,
legend_magnitudes (optional, see `DEFAULT_CATALOG`). Catalogs sized by
magnitude get example markers in a legend beside the globe, placed with
the top-level "legend_position" and "legend_box" keys.

Usage:
    python plot_globe_inset.py globe_inset_example.json --out inset.png
"""
import argparse
import glob
import json
import tempfile
from pathlib import Path

import numpy as np
import pygmt
from obspy import read_events

# Map-level styling, overridden by top-level keys in the JSON file
DEFAULTS = {
    "projection": "orthographic",  # name in PROJECTIONS, or a GMT -J
    "region": None,  # GMT -R; None picks a whole-world region
    "center": None,  # [lon, lat]; None centers on the mean of the points
    "width": "8c",
    "land": "seashell",
    "water": "lightblue",
    "shorelines": "1/0.8p,black",
    "borders": "1/0.25p,black",
    "frame_pen": "1.5p,black",
    "frame": "0",  # GMT -B; "0" is a plain outline, e.g. "afg" annotates
    "transparent": True,
    "parallels": [],
    # GMT legend position and box; default sits right of the globe
    "legend_position": "JMR+jML+o0.3c/0c",
    "legend_box": "+gwhite+p0.5p,black",
}

# Named projections: GMT -J template and the latitude limit of the default
# region (None for -Rg, the whole globe, used by azimuthal projections).
# Azimuthal projections are centered on "center"; the rest only use its
# longitude as the central meridian.
PROJECTIONS = {
    "orthographic": ("G{lon}/{lat}/{width}", None),
    "azimuthal_equidistant": ("E{lon}/{lat}/{width}", None),
    # Horizon of 180 shows the whole world (GMT's default is a hemisphere)
    "lambert_azimuthal": ("A{lon}/{lat}/180/{width}", None),
    "robinson": ("N{lon}/{width}", 90),
    "mollweide": ("W{lon}/{width}", 90),
    "hammer": ("H{lon}/{width}", 90),
    "winkel_tripel": ("R{lon}/{width}", 90),
    "eckert4": ("Kf{lon}/{width}", 90),
    "eckert6": ("Ks{lon}/{width}", 90),
    "equidistant_cylindrical": ("Q{lon}/{width}", 90),
    "miller": ("J{lon}/{width}", 90),
    "mercator": ("M{lon}/{width}", 80),  # poles are at infinity
}

# Per-point styling, overridden by keys on each point in the JSON file
DEFAULT_POINT = {
    "style": "c0.2c",
    "fill": "red",
    "pen": "black",
    "label": None,
    "font": "10p,Helvetica-Bold,black",
    "justify": "LM",  # Left-Middle: text starts to the right of the point
    "offset": "0.25c/0c",  # nudge the text away from the marker
    "arrow": False,  # draw an arrow from the label to the marker
    "arrow_pen": "0.75p,black",
    "arrow_fill": "black",  # arrowhead fill
    "arrow_head": "0.2c",  # arrowhead length
    "arrow_gap": 0.05,  # cm left clear at the label and marker ends
}

# Per-catalog styling, overridden by keys on each catalog in the JSON file
DEFAULT_CATALOG = {
    "symbol": "c",  # GMT symbol code without a size, e.g. c, a, t, s
    "size": 0.1,  # marker size in cm when not sized by magnitude
    "size_by_magnitude": True,
    "size_range": [0.05, 0.3],  # marker size in cm at mag_range ends
    "mag_range": None,  # [min, max]; None uses the catalog's range
    "color_by": None,  # "depth", "time", or None for a single fill
    "fill": "red",  # used when color_by is None
    "cmap": "viridis",
    "cmap_reverse": False,
    "color_range": None,  # [min, max]; None uses the catalog's range
    "colorbar": True,
    "colorbar_label": None,  # None picks "Depth (km)" or "Year"
    "pen": "0.25p,black",
    "name": None,  # legend header for this catalog
    "legend": True,  # show magnitude sizes in the legend
    "legend_magnitudes": None,  # e.g. [4, 5, 6]; None picks round values
}


def load_config(path):
    """
    Read the JSON config and fill in defaults for map, point and catalog
    styling.

    :type path: str
    :param path: path to the JSON file
    :rtype: dict
    :return: config with every key in `DEFAULTS` present, each entry in
        "points" containing every key in `DEFAULT_POINT`, and each entry
        in "catalogs" containing every key in `DEFAULT_CATALOG` plus its
        file paths resolved relative to the JSON file
    """
    with open(path) as f:
        user = json.load(f)

    config = {**DEFAULTS, **user}
    points = []
    for i, point in enumerate(user.get("points", [])):
        if "lat" not in point or "lon" not in point:
            raise ValueError(f"point {i} in {path} is missing 'lat' or 'lon'")
        points.append({**DEFAULT_POINT, **point})
    config["points"] = points

    json_dir = Path(path).parent
    catalogs = []
    for i, catalog in enumerate(user.get("catalogs", [])):
        if "files" not in catalog:
            raise ValueError(f"catalog {i} in {path} is missing 'files'")
        catalog = {**DEFAULT_CATALOG, **catalog}
        if catalog["color_by"] not in (None, "depth", "time"):
            raise ValueError(f"catalog {i} in {path}: color_by must be "
                             "'depth', 'time' or null")
        patterns = catalog["files"]
        if isinstance(patterns, str):
            patterns = [patterns]
        catalog["files"] = []
        for pattern in patterns:
            matches = sorted(glob.glob(str(json_dir / pattern)))
            if not matches:
                raise FileNotFoundError(
                    f"catalog {i} in {path}: no files match {pattern}")
            catalog["files"].extend(matches)
        catalogs.append(catalog)
    config["catalogs"] = catalogs
    return config


def read_catalog(files):
    """
    Read QuakeML files with ObsPy and pull out the preferred (or first)
    origin and magnitude of each event. Events without an origin are
    skipped.

    :type files: list of str
    :param files: QuakeML file paths
    :rtype: dict
    :return: numpy arrays keyed by "lon", "lat", "depth" (km), "time"
        (decimal year) and "mag" (NaN where an event has no magnitude)
    """
    rows = []
    for fid in files:
        for event in read_events(fid):
            origin = event.preferred_origin() or (
                event.origins[0] if event.origins else None)
            if origin is None or origin.latitude is None:
                continue
            magnitude = event.preferred_magnitude() or (
                event.magnitudes[0] if event.magnitudes else None)
            t = origin.time
            year_start = t.__class__(t.year, 1, 1)
            year_len = t.__class__(t.year + 1, 1, 1) - year_start
            rows.append((
                origin.longitude, origin.latitude,
                (origin.depth or 0) / 1000,  # QuakeML depth is in meters
                t.year + (t - year_start) / year_len,
                magnitude.mag if magnitude else np.nan,
            ))
    print(f"Read {len(rows)} events from {len(files)} QuakeML file(s)")
    lon, lat, depth, time, mag = (np.array(col, dtype=float)
                                  for col in zip(*rows)) if rows else \
        [np.array([])] * 5
    return {"lon": lon, "lat": lat, "depth": depth, "time": time,
            "mag": mag}


def magnitude_sizes(mag, size_range, mag_range=None):
    """
    Linearly map magnitudes onto marker sizes, clipping magnitudes
    outside `mag_range`. Events with no magnitude get the smallest size.

    :type mag: numpy.ndarray
    :param mag: event magnitudes, NaN where unknown
    :type size_range: list of float
    :param size_range: [smallest, largest] marker size in cm
    :type mag_range: list of float or None
    :param mag_range: [min, max] magnitude mapped to `size_range`; None
        uses the min and max of `mag`
    :rtype: numpy.ndarray
    :return: marker sizes in cm
    """
    mag_range = resolve_mag_range(mag, mag_range)
    if mag_range is None:
        return np.full_like(mag, size_range[0])
    span = mag_range[1] - mag_range[0]
    scaled = (mag - mag_range[0]) / span if span > 0 else np.ones_like(mag)
    scaled = np.nan_to_num(np.clip(scaled, 0, 1), nan=0.)
    return size_range[0] + scaled * (size_range[1] - size_range[0])


def resolve_mag_range(mag, mag_range=None):
    """
    Pick the magnitude range that marker sizes are scaled over.

    :type mag: numpy.ndarray
    :param mag: event magnitudes, NaN where unknown
    :type mag_range: list of float or None
    :param mag_range: user-set [min, max]; returned unchanged if given
    :rtype: list of float or None
    :return: [min, max] magnitude, or None if no event has a magnitude
    """
    if mag_range is not None:
        return mag_range
    if np.isnan(mag).all():
        return None
    return [float(np.nanmin(mag)), float(np.nanmax(mag))]


def legend_magnitudes(mag_range, max_entries=5):
    """
    Choose round magnitudes spanning `mag_range` to show in the legend:
    whole magnitudes if at least two fall inside the range, otherwise
    the range's ends rounded to one decimal.

    :type mag_range: list of float
    :param mag_range: [min, max] magnitude
    :type max_entries: int
    :param max_entries: most magnitudes to return; whole magnitudes are
        thinned to every 2nd, 3rd, ... value to fit
    :rtype: list of float
    :return: magnitudes in increasing order
    """
    mags = np.arange(np.ceil(mag_range[0]), np.floor(mag_range[1]) + 1)
    if mags.size < 2:
        return sorted({round(mag_range[0], 1), round(mag_range[1], 1)})
    step = int(np.ceil(mags.size / max_entries))
    return list(mags[::-1][::step][::-1])


def magnitude_legend_lines(events, catalog):
    """
    Build GMT legend-spec lines showing marker size for a few
    magnitudes of one catalog.

    :type events: dict
    :param events: output of `read_catalog`
    :type catalog: dict
    :param catalog: catalog styling, see `DEFAULT_CATALOG`
    :rtype: list of str
    :return: legend-spec lines, empty if the catalog isn't sized by
        magnitude, has its legend turned off, or has no magnitudes
    """
    if not (catalog["size_by_magnitude"] and catalog["legend"]):
        return []
    mag_range = resolve_mag_range(events["mag"], catalog["mag_range"])
    if mag_range is None:
        return []
    mags = catalog["legend_magnitudes"] or legend_magnitudes(mag_range)
    sizes = magnitude_sizes(np.array(mags, dtype=float),
                            catalog["size_range"], mag_range)
    # Symbols are drawn in a column wide enough for the largest marker
    x_symbol = max(catalog["size_range"]) / 2 + 0.1
    x_text = max(catalog["size_range"]) + 0.3
    fill = catalog["fill"] if catalog["color_by"] is None else "gray70"

    lines = []
    if catalog["name"]:
        lines.append(f"H 9p,Helvetica-Bold {catalog['name']}")
    for mag, size in zip(mags, sizes):
        lines.append(f"S {x_symbol:g}c {catalog['symbol']} {size:.3f}c "
                     f"{fill} {catalog['pen']} {x_text:g}c M{mag:g}")
    return lines


def plot_catalog(fig, events, catalog, width_cm, n_colorbars):
    """
    Plot one catalog's events, sized by magnitude and colored by depth or
    time if requested, with a colorbar below the globe.

    :type fig: pygmt.Figure
    :param fig: figure to draw on
    :type events: dict
    :param events: output of `read_catalog`
    :type catalog: dict
    :param catalog: catalog styling, see `DEFAULT_CATALOG`
    :type width_cm: float
    :param width_cm: globe width in cm, used to size the colorbar
    :type n_colorbars: int
    :param n_colorbars: colorbars already drawn, used to stack this one
        beneath them
    :rtype: bool
    :return: True if a colorbar was drawn
    """
    if events["lon"].size == 0:
        return False

    if catalog["size_by_magnitude"]:
        sizes = magnitude_sizes(events["mag"], catalog["size_range"],
                                catalog["mag_range"])
    else:
        sizes = np.full_like(events["lon"], catalog["size"])
    # Draw the biggest events first so smaller ones aren't hidden
    order = np.argsort(-sizes)
    style = f"{catalog['symbol']}c"  # sizes come from `size`, in cm

    color_by = catalog["color_by"]
    if color_by is None:
        fig.plot(x=events["lon"][order], y=events["lat"][order],
                 size=sizes[order], style=style, fill=catalog["fill"],
                 pen=catalog["pen"])
        return False

    values = events[color_by]
    vmin, vmax = catalog["color_range"] or [values.min(), values.max()]
    if vmin == vmax:
        vmin, vmax = vmin - 0.5, vmax + 0.5
    pygmt.makecpt(cmap=catalog["cmap"], series=[vmin, vmax],
                  reverse=catalog["cmap_reverse"])
    fig.plot(x=events["lon"][order], y=events["lat"][order],
             size=sizes[order], fill=values[order], cmap=True,
             style=style, pen=catalog["pen"])
    if catalog["colorbar"]:
        label = catalog["colorbar_label"] or {
            "depth": "Depth (km)", "time": "Year"}[color_by]
        offset = 0.5 + n_colorbars * 1.3
        fig.colorbar(position=(f"JBC+w{0.8 * width_cm:g}c/0.2c"
                               f"+o0c/{offset:g}c+h"),
                     frame=[f"x+l{label}"])
        return True
    return False


def mean_center(lats, lons):
    """
    Find the geographic mean of a set of locations, used as the globe's
    projection center. Averages unit vectors so longitudes that straddle
    the antimeridian (e.g. 179 and -179) average correctly.

    :type lats: array-like
    :param lats: latitudes in degrees
    :type lons: array-like
    :param lons: longitudes in degrees
    :rtype: tuple of float
    :return: (lon, lat) in degrees
    """
    lats = np.radians(lats)
    lons = np.radians(lons)
    x = np.mean(np.cos(lats) * np.cos(lons))
    y = np.mean(np.cos(lats) * np.sin(lons))
    z = np.mean(np.sin(lats))
    lon = np.degrees(np.arctan2(y, x))
    lat = np.degrees(np.arctan2(z, np.hypot(x, y)))
    return lon, lat


def plot_globe(config):
    """
    Draw the global map with parallels, catalog events and points from
    `config`.

    :type config: dict
    :param config: output of `load_config`
    :rtype: pygmt.Figure
    :return: the finished figure
    """
    catalog_events = [read_catalog(c["files"]) for c in config["catalogs"]]

    lats = [p["lat"] for p in config["points"]]
    lons = [p["lon"] for p in config["points"]]
    for events in catalog_events:
        lats.extend(events["lat"])
        lons.extend(events["lon"])
    if config["center"] is not None:
        center_lon, center_lat = config["center"]
    elif lats:
        center_lon, center_lat = mean_center(lats, lons)
    else:
        center_lon, center_lat = 0, 0

    projection, region = build_projection(config, center_lon, center_lat)
    frame_size = map_size(projection, region)

    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_PEN=config["frame_pen"])
    fig.coast(
        projection=projection,
        region=region,
        land=config["land"],
        water=config["water"],
        shorelines=config["shorelines"],
        borders=config["borders"],
    )
    # Draw the border/frame around the map
    with pygmt.config(MAP_FRAME_TYPE="plain"):
        fig.basemap(frame=config["frame"])

    # Lines of constant latitude, e.g. the Arctic Circle
    circle_lons = np.arange(-180, 180.1, 1)
    for parallel in config["parallels"]:
        fig.plot(x=circle_lons,
                 y=np.full_like(circle_lons, parallel["lat"]),
                 pen=parallel.get("pen", "1p,black,--"))

    # Catalogs go beneath the points so labeled sites stay visible
    n_colorbars = 0
    for events, catalog in zip(catalog_events, config["catalogs"]):
        n_colorbars += plot_catalog(fig, events, catalog, frame_size[0],
                                    n_colorbars)

    legend_lines = []
    for events, catalog in zip(catalog_events, config["catalogs"]):
        lines = magnitude_legend_lines(events, catalog)
        if lines and legend_lines:
            legend_lines.append("G 0.15c")  # gap between catalogs
        legend_lines.extend(lines)
    if legend_lines:
        with tempfile.NamedTemporaryFile("w", suffix=".txt") as spec:
            spec.write("\n".join(legend_lines) + "\n")
            spec.flush()
            fig.legend(spec=spec.name, position=config["legend_position"],
                       box=config["legend_box"] or False)

    for point in config["points"]:
        fig.plot(x=point["lon"], y=point["lat"], style=point["style"],
                 fill=point["fill"], pen=point["pen"])
        if point["label"]:
            draw_label(fig, point)

    # Arrows are drawn last, on a Cartesian frame (in cm) laid exactly over
    # the map, because the label offset is a distance on the page
    for point in config["points"]:
        if point["label"] and point["arrow"]:
            draw_label_arrow(fig, point, projection, region, frame_size)
    return fig


def draw_label(fig, point, line_spacing=1.2):
    """
    Draw a point's label. GMT text doesn't break lines on "\\n", so a
    multi-line label is split and each line drawn with its own vertical
    offset, stacked so the block as a whole honors the vertical part of
    `justify` (T: hangs down from the anchor, M: centered on it, B: sits
    on it).

    :type fig: pygmt.Figure
    :param fig: figure to draw on
    :type point: dict
    :param point: point styling, see `DEFAULT_POINT`
    :type line_spacing: float
    :param line_spacing: line height as a multiple of the font size
    """
    lines = str(point["label"]).split("\n")
    if len(lines) == 1:
        fig.text(x=point["lon"], y=point["lat"], text=lines[0],
                 font=point["font"], justify=point["justify"],
                 offset=point["offset"])
        return

    dx, dy = (to_cm(v) for v in str(point["offset"]).split("/"))
    font_size = point["font"].split(",")[0]
    height = to_cm(font_size if font_size[-1].isalpha() else
                   f"{font_size}p") * line_spacing
    vertical = next((c for c in point["justify"].upper() if c in "TMB"),
                    "M")
    n = len(lines)
    for i, line in enumerate(lines):
        if vertical == "T":
            shift = -i * height
        elif vertical == "B":
            shift = (n - 1 - i) * height
        else:
            shift = ((n - 1) / 2 - i) * height
        fig.text(x=point["lon"], y=point["lat"], text=line,
                 font=point["font"], justify=point["justify"],
                 offset=f"{dx:g}c/{dy + shift:g}c")


def to_cm(value):
    """
    Convert a GMT length such as "0.25c", "0.1i" or "6p" to centimeters.
    Plain numbers are taken to be centimeters.

    :type value: str or float
    :param value: GMT length
    :rtype: float
    :return: length in cm
    """
    value = str(value).strip()
    units = {"c": 1, "i": 2.54, "p": 2.54 / 72}
    if value and value[-1] in units:
        return float(value[:-1]) * units[value[-1]]
    return float(value)


def build_projection(config, center_lon, center_lat):
    """
    Turn the config's "projection" into a GMT -J string and matching -R
    region. Names in `PROJECTIONS` are filled in with the map center and
    width; anything else is passed to GMT unchanged as a full -J string
    (e.g. "G-150/45/10c").

    :type config: dict
    :param config: output of `load_config`
    :type center_lon: float
    :param center_lon: projection center longitude in degrees
    :type center_lat: float
    :param center_lat: projection center latitude in degrees (only used
        by azimuthal projections)
    :rtype: tuple of str
    :return: (projection, region)
    """
    name = config["projection"]
    if name in PROJECTIONS:
        template, max_lat = PROJECTIONS[name]
        projection = template.format(lon=f"{center_lon:g}",
                                     lat=f"{center_lat:g}",
                                     width=config["width"])
    else:
        projection, max_lat = name, None
    if config["region"]:
        region = config["region"]
    elif max_lat is None:
        region = "g"
    else:
        # Whole-world map with the center longitude in the middle
        region = (f"{center_lon - 180:g}/{center_lon + 180:g}/"
                  f"{-max_lat:g}/{max_lat:g}")
    return projection, region


def run_mapproject(args):
    """
    Run GMT mapproject and return its numeric output.

    :type args: str
    :param args: mapproject arguments, without the output redirect
    :rtype: numpy.ndarray
    :return: 2D array of the output rows (empty if there were none)
    """
    with tempfile.NamedTemporaryFile(suffix=".txt") as out:
        with pygmt.clib.Session() as lib:
            lib.call_module("mapproject", f"{args} ->{out.name}")
        return np.loadtxt(out.name, ndmin=2)


def map_size(projection, region):
    """
    Get the width and height of the map frame on the page.

    :type projection: str
    :param projection: GMT -J string
    :type region: str
    :param region: GMT -R string
    :rtype: tuple of float
    :return: (width, height) in cm
    """
    width, height = run_mapproject(f"-J{projection} -R{region} -W")[0]
    return width, height


def project_xy(lon, lat, projection, region):
    """
    Project a location onto the page, measured from the lower-left
    corner of the map frame.

    :type lon: float
    :param lon: longitude in degrees
    :type lat: float
    :param lat: latitude in degrees
    :type projection: str
    :param projection: GMT -J string
    :type region: str
    :param region: GMT -R string
    :rtype: tuple
    :return: (x, y) in cm, or None if the location isn't on the map (e.g.
        the far side of an orthographic globe)
    """
    with tempfile.NamedTemporaryFile("w", suffix=".txt") as src:
        src.write(f"{lon} {lat}\n")
        src.flush()
        # -S drops locations outside the map
        xy = run_mapproject(f"{src.name} -J{projection} -R{region} -S")
    return tuple(xy[0]) if xy.size else None


def draw_label_arrow(fig, point, projection, region, frame_size):
    """
    Draw an arrow from a point's offset label to its marker. The arrow
    starts at the label's justification anchor (e.g. the left-middle of
    the text for justify="LM") and stops short of the marker's edge.

    :type fig: pygmt.Figure
    :param fig: figure to draw on
    :type point: dict
    :param point: point styling, see `DEFAULT_POINT`
    :type projection: str
    :param projection: GMT -J string of the map
    :type region: str
    :param region: GMT -R string of the map
    :type frame_size: tuple of float
    :param frame_size: (width, height) of the map frame in cm
    """
    marker = project_xy(point["lon"], point["lat"], projection, region)
    if marker is None:
        print(f"skipping arrow for {point['label']}: point is not visible "
              "on the map")
        return
    dx, dy = (to_cm(v) for v in str(point["offset"]).split("/"))
    length = np.hypot(dx, dy)
    # Marker size is the number after the symbol code, e.g. "c0.2c"
    marker_radius = to_cm(point["style"][1:]) / 2
    start_gap = point["arrow_gap"]
    end_gap = marker_radius + point["arrow_gap"]
    if length <= start_gap + end_gap:
        print(f"skipping arrow for {point['label']}: offset is too small "
              "to fit an arrow between the label and marker")
        return

    ux, uy = -dx / length, -dy / length  # unit vector, label -> marker
    x0 = marker[0] + dx + ux * start_gap
    y0 = marker[1] + dy + uy * start_gap
    x1 = marker[0] + ux * -end_gap
    y1 = marker[1] + uy * -end_gap
    fig.plot(x=[x0], y=[y0], direction=([x1], [y1]),
             style=f"v{point['arrow_head']}+e+s+h0.5",
             pen=point["arrow_pen"], fill=point["arrow_fill"],
             projection=f"X{frame_size[0]:g}c/{frame_size[1]:g}c",
             region=[0, frame_size[0], 0, frame_size[1]])


def parse_args():
    """
    Parse command line arguments.

    :rtype: argparse.Namespace
    :return: parsed arguments
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", help="JSON file of points and styling")
    parser.add_argument("--out", default="globalinset.png",
                        help="output figure path (default: globalinset.png)")
    parser.add_argument("--show", action="store_true",
                        help="open the figure after saving")
    return parser.parse_args()


def main():
    """Load the config, draw the globe, and save it."""
    args = parse_args()
    config = load_config(args.config)
    fig = plot_globe(config)
    fig.savefig(args.out, transparent=config["transparent"])
    print(f"Saved figure to {args.out}")
    if args.show:
        fig.show()


if __name__ == "__main__":
    main()
