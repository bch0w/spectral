"""Create a PyGMT map for a bounding box and add markers for stations/cities.

This script plots topography, contour elevations, coastlines, and allows
simple latitude/longitude markers with optional labels.
"""
import os
from typing import List, Dict, Optional

import pygmt


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
) -> pygmt.Figure:
    """Plot a bounding-box map with topography, contours, coastlines, and markers.

    Parameters
    ----------
    region
        [minlon, maxlon, minlat, maxlat]
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
        - lon: longitude
        - lat: latitude
        - style: GMT symbol style string (default "c0.4c")
        - color: fill color
        - pen: outline pen
        - label: optional text label
    title
        Optional title text.
    show
        If True, show the figure after saving.

    Returns
    -------
    pygmt.Figure
        The generated figure.
    """
    fig = pygmt.Figure()

    # Load satellite-derived topography/bathymetry from the GMT remote repository.
    # Use the GEBCO-derived dataset for higher-quality relief imagery.
    topo_grid = pygmt.datasets.load_earth_relief(
        resolution=topo_resolution,
        region=region,
        data_source="gebco",
    )

    fig.grdimage(
        grid=topo_grid,
        projection=projection,
        region=region,
        cmap="gmt/geo",
        shading=True,
        frame=["a", "+t" + title if title else ""],
    )

    # Contour the topography
    fig.grdcontour(
        grid=topo_grid,
        interval=contour_interval,
        pen="0.5p,black,gray",
        annotation=f"{contour_interval}+f8p",
        limit=[-6000, 9000],
    )

    # Add coastlines over the topo layer.
    fig.coast(
        region=region,
        projection=projection,
        shorelines=["1/0.8p,black"],
        frame=["a", "+t" + title if title else ""],
        resolution="i",
    )

    # Plot markers and optional labels
    if markers:
        for marker in markers:
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


if __name__ == "__main__":
    # Example usage: bounding box around Wellington, NZ
    region_example = [125.5, 130.5, 34, 42.5]
    markers_example = [
        {"lon": 129.0297, "lat": 41.33, "label": "NKNTS", "color": "red"},
        # {"lon": 175.28, "lat": -40.90, "label": "Upper Hutt", "color": "green"},
        # {"lon": 175.83, "lat": -41.28, "label": "Kapiti", "color": "magenta"},
    ]

    plot_bbox_map(
        region=region_example,
        output="mapping/pygmt/bbox_map_example.png",
        topo_resolution="01m",
        projection="M12c",
        contour_interval=200,
        markers=markers_example,
        title="Test",
    )
