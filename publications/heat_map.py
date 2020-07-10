"""
Make a heat map of event locations
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from obspy import read_events
from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import beach
from pyatoa.visuals.map_maker import event_beachball
from pyatoa.utils.calculate import normalize_a_to_b, myround

from big_map import (place_scalebar, initiate_basemap, draw_domain_bounds)


def get_lat_lons(fid):
    """
    Get latitude and longitude values from a catalog
    """
    cat = read_events(fid)

    # First get some array info for scale bars and colormap
    print(f"{len(cat)} events in catalog")
    lats, lons = [], []
    for i, event in enumerate(cat):
        try:
            origin = event.preferred_origin()
            mag = event.preferred_magnitude().mag
        except AttributeError:
            origin = event.origins[0]
            mag = event.magnitudes[0].mag

        if (origin.depth * 1E-3 > 60) or (mag < 4.4):
            continue

        lats.append(origin.latitude)
        lons.append(origin.longitude)

    return lats, lons


def draw_heatmap(m, nx, ny, lats, lons, map_corners, cmap="Oranges"):
    """
    Draw a heatmap on top of a basemap
    http://qingkaikong.blogspot.com/2016/02/
                                plot-earthquake-heatmap-on-basemap-and.html

    :type nx: int
    :param nx: number of bins for the x-axis
    :type ny: int
    :param ny: number of bins for the y-axis
    :type lats: np.array
    :param lats: array of latitude values
    :type lons: array of longitude values
    :type map_corners: dict
    :param map_corners: dict of map latitude and longitude values
    """
    # form the bins
    lon_bins = np.linspace(map_corners["lon_min"], map_corners["lon_max"], nx)
    lat_bins = np.linspace(map_corners["lat_min"], map_corners["lat_max"], ny)
        
    # aggregate the number of earthquakes in each bin, we will only use the density
    density, lat_edges, lon_edges = np.histogram2d(lats, lons, 
                                                   [lat_bins, lon_bins])

    # get the mesh for the lat and lon
    lon_bins_2d, lat_bins_2d = np.meshgrid(lon_bins, lat_bins)

    # convert the bin mesh to map coordinates:
    xs, ys = m(lon_bins_2d, lat_bins_2d) # will be plotted using pcolormesh

    # Here adding one row and column at the end of the matrix, so that 
    # density has same dimension as xs, ys, otherwise, using shading='gouraud'
    # will raise error
    density = np.hstack((density,np.zeros((density.shape[0], 1))))
    density = np.vstack((density,np.zeros((density.shape[1]))))

    # Plot heatmap with the custom color map
    pcm = plt.pcolormesh(xs, ys, density, cmap=cmap, shading="flat")
    pcm.set_clim(0, 5.)

    # Add color bar and 
    cbar = plt.colorbar(extend='max', shrink=0.2, aspect=10)
    cbar.set_label('Number of earthquakes', size=18)

    x, y = m(lons, lats)
    plt.scatter(x, y, marker="o", c="r", edgecolor="r", linewidth=1.)

    # For plotting the grid points
    # plt.scatter(xs, ys, marker="s", c="None", edgecolor="k")


if __name__ == "__main__":
    cat_fid = ("/Users/Chow/Documents/academic/vuw/data/events/"
               "nz_north_cat_253.xml")

    # Figure kwargs
    figsize = (10, 12)
    dpi = 100

    # Map kwargs
    map_corners={'lat_min': -42.75, 'lat_max': -36.75,  
                 'lon_min': 172.5, 'lon_max': 178.8}
    domain_bounds={'lat_min': -42.5, 'lat_max': -37.0,  
                   'lon_min': 173.0, 'lon_max': 178.5}
    coastline_linewidth = 2
    continent_color= "None"
    lake_color = "None"
    fill_color = "None"
    scalebar_location = "lower-right"
    scalebar_fontsize = 18.
    scalebar_linewidth = 2.75
    fontsize = 22.5

    # Heat map kwargs
    nx = 10
    ny = 12
    cmap = "Greys"

    # Plotting starts here
    f = plt.figure(figsize=figsize, dpi=dpi)
    ax = f.add_subplot(111)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(2)

    m = initiate_basemap(map_corners, continent_color=continent_color,
                         lake_color=lake_color, fill_color=fill_color,
                         coastline_linewidth=coastline_linewidth, scalebar=True,
                         scalebar_location=scalebar_location, fontsize=fontsize,
                         scalebar_fontsize=scalebar_fontsize, 
                         scalebar_linewidth=scalebar_linewidth
                         )


    draw_domain_bounds(m, domain_bounds)
    lats, lons = get_lat_lons(cat_fid)
    draw_heatmap(m, nx, ny, lats, lons, map_corners, cmap)

    # make image bigger:
    # plt.gcf().set_size_inches(12,12)

    f.tight_layout()
    plt.savefig("./figures/heat_map.png")
    plt.show()
