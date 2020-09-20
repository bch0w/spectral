"""
Make a global projection map using basemap, for reference to NZ geographic loc
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def initiate_basemap(lat_0, lon_0, **kwargs):
    """
    set up the basemap object in the same way each time

    :type map_center: dict of floats
    :param map_center: {lat_0, lon_0}
    :rtype m: Basemap
    :return m: basemap object
    """
    continent_color = kwargs.get("continent_color", "w")
    lake_color = kwargs.get("lake_color", "w")
    coastline_zorder = kwargs.get("coastline_zorder", 5)
    coastline_linewidth = kwargs.get("coastline_linewidth", 1.0)
    fill_color = kwargs.get("fill_color", "w")
    fontsize = kwargs.get("fontsize", 13)
    resolution = kwargs.get("resolution", "h")

    # Initiate map and draw in style
    m = Basemap(projection="ortho", resolution=resolution, lat_0=lat_0, 
                lon_0=lon_0)
    m.drawcoastlines(linewidth=coastline_linewidth, zorder=coastline_zorder)
    m.fillcontinents(color=continent_color, lake_color=lake_color)
    m.drawmapboundary(fill_color=fill_color)
    m.drawparallels(np.arange(-90., 120., 30.))
    m.drawmeridians(np.arange(0., 420., 60.))

    return m

def draw_domain_bounds(m, bounds):
    """
    Outline a domain bound on the basemap
    """
    lats = [bounds["lat_min"], bounds["lat_max"],
            bounds["lat_max"], bounds["lat_min"]]
    lons = [bounds["lon_min"], bounds["lon_min"],
            bounds["lon_max"], bounds["lon_max"]]
    x, y = m(lons, lats)
    xy = zip(x,y)
    poly = mpl.patches.Polygon(list(xy), facecolor='None', edgecolor="r",
                               linewidth=1.5, linestyle='-', alpha=1)

    plt.gca().add_patch(poly)


if __name__ == "__main__":
    # Cook Strait: -41.3912, 174.5222
    m = initiate_basemap(lat_0=-41.3912, lon_0=174.5222, resolution="l")
    draw_domain_bounds(m, bounds={"lat_min": -42.5, "lat_max": -37.0,
                                  "lon_min": 173.0, "lon_max": 178.5}
                       )
    plt.show()

