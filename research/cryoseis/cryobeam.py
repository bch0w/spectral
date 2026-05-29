"""
Beamformer tuned for 2024 Gulkanaseis deployment
Beamforming source code inspired by: https://github.com/schipp/fast_beamforming
"""
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Transformer
from typing import List, Dict, Optional, Tuple


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


def array_geometry(origin_latlon=None, plot=True):
    """
    Define array geometry spatially and plot

    Parameters
    ----------
    origin (tuple of float)
        Arbitrarily defined origin (0, 0) to ensure UTM coordinates are not 
        very large values. By default will set to middle station location.
        Alterantively I have chosen one arbitrary center to be the rock face 
        east of Gabriel Ice Fall (63.272544, -145.431690)
    """
    # Gulkana Seismic Deployment 2024
    stations = {
        "100": (63.265612, -145.415025),
        "101": (63.265362, -145.41407),
        "102": (63.265448, -145.41381),
        "103": (63.265203, -145.415355),
        "104": (63.265062, -145.414323),
        "105": (63.26511,	-145.413887),
        "106": (63.26495,	-145.415472),
        "107": (63.264683, -145.414657),
        "108": (63.264878, -145.414107),
        "109": (63.26392,	-145.41601),
        "110": (63.26466,	-145.41557),
        "111": (63.264278, -145.415858),
    }    

    if origin_latlon is None:
        origin_latlon = stations["107"]

    origin_lat, origin_lon = origin_latlon
    origin_x, origin_y = lonlat_to_utm(lon=origin_lon, lat=origin_lat)[:2]
    
    # Convert to UTM coordinates so we are in a Cartesian reference frame
    stations_utm  =  {}
    for name, coords in stations.items():
        lat, lon = coords
        x, y = lonlat_to_utm(lon, lat)[:2]
        x -= origin_x
        y -= origin_y
        stations_utm[name] = (x, y)

        if plot:
            plt.scatter(x, y, c="k")
            plt.text(x, y, name)

    if plot:
        plt.xlabel("Easting [m]")
        plt.ylabel("Northing [m]")
        plt.title("Gulkanaseis 2024")
        plt.show()

    return stations_utm



def main():
    array_geometry()

if __name__ == "__main__":
    main()