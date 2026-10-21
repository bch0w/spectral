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
offset (optional, see `DEFAULT_POINT`). If "center" is omitted the globe
is centered on the mean location of the points and catalog events.

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
cmap_reverse, color_range, colorbar, colorbar_label, pen (optional, see
`DEFAULT_CATALOG`).

Usage:
    python plot_globe_inset.py globe_inset_example.json --out inset.png
"""
import argparse
import glob
import json
from pathlib import Path

import numpy as np
import pygmt
from obspy import read_events

# Map-level styling, overridden by top-level keys in the JSON file
DEFAULTS = {
    "center": None,  # [lon, lat]; None centers on the mean of the points
    "width": "8c",
    "land": "seashell",
    "water": "lightblue",
    "shorelines": "1/0.8p,black",
    "borders": "1/0.25p,black",
    "frame_pen": "1.5p,black",
    "transparent": True,
    "parallels": [],
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
    if mag_range is None:
        if np.isnan(mag).all():
            return np.full_like(mag, size_range[0])
        mag_range = [np.nanmin(mag), np.nanmax(mag)]
    span = mag_range[1] - mag_range[0]
    scaled = (mag - mag_range[0]) / span if span > 0 else np.ones_like(mag)
    scaled = np.nan_to_num(np.clip(scaled, 0, 1), nan=0.)
    return size_range[0] + scaled * (size_range[1] - size_range[0])


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
    Draw the orthographic globe with parallels, catalog events and points
    from `config`.

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

    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_PEN=config["frame_pen"])
    fig.coast(
        # Orthographic projection (G) centered on the region of interest
        projection=f"G{center_lon:g}/{center_lat:g}/{config['width']}",
        region="g",  # global
        land=config["land"],
        water=config["water"],
        shorelines=config["shorelines"],
        borders=config["borders"],
    )
    # Draw the border/frame around the globe
    fig.basemap(frame=True)

    # Lines of constant latitude, e.g. the Arctic Circle
    circle_lons = np.arange(-180, 180.1, 1)
    for parallel in config["parallels"]:
        fig.plot(x=circle_lons,
                 y=np.full_like(circle_lons, parallel["lat"]),
                 pen=parallel.get("pen", "1p,black,--"))

    # Catalogs go beneath the points so labeled sites stay visible
    width_cm = float(str(config["width"]).rstrip("c"))
    n_colorbars = 0
    for events, catalog in zip(catalog_events, config["catalogs"]):
        n_colorbars += plot_catalog(fig, events, catalog, width_cm,
                                    n_colorbars)

    for point in config["points"]:
        fig.plot(x=point["lon"], y=point["lat"], style=point["style"],
                 fill=point["fill"], pen=point["pen"])
        if point["label"]:
            fig.text(x=point["lon"], y=point["lat"], text=point["label"],
                     font=point["font"], justify=point["justify"],
                     offset=point["offset"])
    return fig


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
