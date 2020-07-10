"""
Make a map with events and receivers, plate interface countours and landmarks
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import beach
from obspy import read_inventory, read_events
from pyatoa.utils.calculate import normalize_a_to_b, myround
from pyatoa.visuals.map_maker import event_beachball, interpolate_and_contour

from big_map import (place_scalebar, initiate_basemap, draw_domain_bounds,
                     plot_stations, plot_beachballs)


def plot_cities(m, names=False, **kwargs):
    """
    Plot major population centers
    
    qwer
    
    :type m: basemap
    :param m: map to plot cities on
    :type names: bool
    :param names: annotate names on map
    """
    markersize = kwargs.get("markersize", 80)
    color = kwargs.get("color", "w")
    fontsize = kwargs.get("fontsize", 12)
    zorder = kwargs.get("zorder", 100)
    alpha = kwargs.get("alpha", 0.85)
    
    city = True
    location = False

    locations = [(-41.28664, 174.77557, "Wellington", city),
                  (-39.4928, 176.9120, "Napier", city),
                  (-42.416665, 173.6833306, "Kaikoura", city),
                  (-38.6857, 176.0702, "Taupo", city),
                  (-39.2968, 174.0634, "Mt. Taranaki", city),
                  (-37.8, 176., "Bay of Plenty", location),
#                   (-39.3314, 177.5017, "Hawke Bay", location),
#                   (-38.6623, 178.0176, "Gisborne", city),
#                   (-40.0741, 174.8985, "Taranaki Bight", location),
#                   (-36.8485, 174.7633, "Auckland", city),
                  (-41.6, 174.3121, "Cook Strait", location),
#                   (-41.1383, 174.0870, "Marlborough Sounds", location)
              (-38.8, 173., "Australian Plate", location),
              (-42., 176.5, "Pacific Plate", location)
             ]

    for loc in locations:
        xy = m(loc[1], loc[0])
        if loc[3]:
            m.scatter(xy[0], xy[1], marker='o', s=markersize, edgecolor='k', 
                      c=color, linewidth=2, zorder=zorder)
            plt.text(xy[0] + 0.01 * (m.xmax - m.xmin), 
                     xy[1] + 0.01 * (m.ymax - m.ymin),
                     s=loc[2], fontsize=fontsize, zorder=zorder, 
                     fontweight='normal', bbox=dict(facecolor='w', fill=True, 
                                                    edgecolor='k', 
                                                    linewidth=1.5, alpha=alpha, 
                                                    boxstyle="round")
                     )
        else:
            plt.text(xy[0] + 0.01 * (m.xmax - m.xmin), 
                     xy[1] + 0.01 * (m.ymax - m.ymin),
                     s = loc[2], fontsize=fontsize, zorder=zorder, 
                     fontweight='normal', bbox=dict(facecolor='w', fill=True, 
                                                    edgecolor='k', 
                                                    linewidth=1.5, alpha=alpha, 
                                                    boxstyle="square")
                     )


def plot_interface_contour(m, fid, **kwargs):
    """
    Plot the plate interface as contours, inline labels for depths
    
    :type m: basemap
    :param m: map convert coordinates
    :type plate: np.array
    :param plate: info read in from plate interface npy file
    :type xi, yi, zi: np.array
    :param xi, yi, zi: gridded data for use in contour
    """
    zorder = kwargs.get("zorder", 100)
    levels = kwargs.get("levels", [0,5,10,15,25,50,100,200,300])
    color = kwargs.get("color", "k")
    linestyle = kwargs.get("linestyle", "dotted")
    linewidth = kwargs.get("linewidth", 3)
    alpha = kwargs.get("alpha", 0.5)
    fontsize = kwargs.get("fontsize", 12)
    fontweight = kwargs.get("fontweight", "normal")
    format_ = kwargs.get("format", '%.0f km')

    plate = np.load(fid)  

    lats = plate[:, 0]
    lons = plate[:, 1]
    depths = plate[:, 2]
    depths = [-1 * _ for _ in depths]
    
    x, y = m(lons, lats)
    xi = np.linspace(min(x), max(x), 1000)
    yi = np.linspace(min(y), max(y), 1000)
    
    xi, yi = np.meshgrid(xi, yi)
    zi = griddata((x, y), depths, (xi, yi))

    cs = m.contour(xi, yi, zi, levels, vmin=0, colors=color, 
                   linestyles=linestyle, alpha=alpha, zorder=zorder, 
                   linewidths=linewidth)
    clabels = plt.clabel(cs, levels[1:], fontsize=fontsize, fmt=format_)

    [txt.set_bbox(dict(facecolor='white', edgecolor='None', pad=0)) 
                                                            for txt in clabels]


if __name__ == "__main__":
    cat_fid = ("/Users/Chow/Documents/academic/vuw/data/events/"
               "checkerboard_30.xml")
    sta_fid = ("/Users/Chow/Documents/academic/vuw/data/specfem/stations/"
               "STATIONS_CHKBD")
    plate_fid = ("/Users/Chow/Documents/academic/vuw/data/carto/"
                 "fault_coordinates/hikurangi_plate_interface.npy")

    # Figure kwargs
    figsize = (10, 12)
    dpi = 100

    # Map kwargs
    map_corners={'lat_min': -42.75, 'lat_max': -36.75,  
                 'lon_min': 172.5, 'lon_max': 178.8}
    domain_bounds={'lat_min': -42.5, 'lat_max': -37.0,  
                   'lon_min': 173.0, 'lon_max': 178.5}
    coastline_linewidth = 2
    continent_color= "whitesmoke"
    lake_color = "azure"
    fill_color = "azure"
    scalebar_location = "lower-right"
    scalebar_fontsize = 18.
    scalebar_linewidth = 2.75
    fontsize = 22.5

    # Station kwargs
    station_markersize = 90
    station_zorder = 100
    station_color = "w"
    station_marker = "v"
    station_edgecolor = "k"

    # Event kwargs
    beachball_cmap = plt.cm.jet_r
    beachball_norm_a = 1.25E4
    beachball_norm_b = 2.75E4
    beachball_linewidths = 1.6
    beachball_zorder = station_zorder - 2
    mag_scale = [4.5, 6.0]

    cbar_shrink = 0.2
    cbar_fontsize = 15
    cbar_tickfontsize = 15
    cbar_labelpad = 17.5

    # City kwargs
    city_zorder = station_zorder-1
    city_markersize = 125
    city_fontsize = 14
    city_names = True
    city_color = "yellow"
    city_alpha = 0.5

    # Plate interface kwargs
    contour_zorder = station_zorder - 5
    contour_levels = [0, 5, 10, 15, 25, 50, 100, 200, 300]
    contour_color = "k"
    contour_linestyle = ":"
    contour_alpha = 1.
    contour_linewidth = 1.5
    clabel_fontsize = 15
    clabel_fontweight = "heavy"
    clabel_format = '%.0f km'
    

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

    plot_stations(m, sta_fid, markersize=station_markersize, 
              zorder=station_zorder, color=station_color, 
              station_edgecolor=station_edgecolor,  
              station_marker=station_marker
              )

    plot_beachballs(m, cat_fid, cmap=beachball_cmap, norm_a=beachball_norm_a, 
                    norm_b=beachball_norm_b, linewidth=beachball_linewidths, 
                    zorder=beachball_zorder, cbar_shrink=cbar_shrink,
                    mag_scale=mag_scale
                    )


    plot_interface_contour(m, plate_fid, zorder=contour_zorder, 
                           levels=contour_levels, color=contour_color, 
                           linestyle=contour_linestyle, alpha=contour_alpha, 
                           fontsize=0, format_=clabel_format, 
                           linewidth=contour_linewidth)    

    plot_cities(m, markersize=city_markersize, zorder=city_zorder,
                names=city_names, color=city_color, alpha=city_alpha, 
                fontsize=city_fontsize)

    f.tight_layout()
    plt.savefig("./figures/rich_map.png")
    plt.show()
