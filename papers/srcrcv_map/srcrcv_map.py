"""
Overly-detailed publication quality source receiver map
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from obspy import read_inventory, read_events

from pyatoa.utils.tools.calculate import normalize_a_to_b, myround
from pyatoa.utils.visuals.map_tools import initiate_basemap, event_beachball, \
        interpolate_and_contour
from pyatoa.utils.visuals.map_plugins import plot_geonet_active_faults

mpl.rcParams['axes.linewidth'] = 2.


def plot_inv(m, inv, **kwargs):
    """
    plot stations from an obspy inventory object
    """
    # place changeable variables up high for quick change
    markersize = kwargs.get("markersize", 100)
    zorder = kwargs.get("zorder", 100)

    x, y, duration = [], [], []
    for net in inv:
        for sta in net:
            name = f"{net.code}.{sta.code}"

            # Determine station up-time
            if sta.end_date is None:
                duration.append(0)
            else:
                duration.append(sta.end_date - sta.start_date)
            xy = m(sta.longitude, sta.latitude)
            x.append(xy[0])
            y.append(xy[1])

    # replace the permanent stations with some value
    for i, dur in enumerate(duration):
        print(dur)
        if not dur:
            duration[i] = max(duration) * 2

    m.scatter(x, y, marker=11, s=markersize, edgecolor='k', c=duration,
              linestyle='-', linewidth=1.25, cmap='viridis', 
              zorder=zorder)


def plot_stations(m, stations_fid, **kwargs):
    """
    plot stations from a Specfem STATION file
    """
    # place changeable variables up high for quick change
    markersize = kwargs.get("markersize", 100)
    zorder = kwargs.get("zorder", 100)

    station_info = np.loadtxt(stations_fid, usecols=(0,1,2,3), dtype=str)
    x, y, colors = [], [], []
    for station in station_info:
        sta, net, lat, lon = station
        name = f"{net}.{sta}"
        xy = m(float(lon), float(lat))
        x.append(xy[0])
        y.append(xy[1])
        if net == "NZ":
            colors.append("k")
        else:
            colors.append("gray")

    m.scatter(x, y, marker=11, s=markersize, edgecolors='k', 
              c=colors, linestyle='-', linewidth=2., zorder=zorder)


def plot_cities(m, names=False, **kwargs):
    """
    Plot major population centers
    """
    markersize = kwargs.get("markersize", 80)
    fontsize = kwargs.get("fontsize", 12)
    zorder = kwargs.get("zorder", 100)

    cities = [(-41.28664, 174.77557, "Wellington"),
              (-39.6381, 176.84918, "Hastings"),
              (-42.416665, 173.6833306, "Kaikoura")]
    
    for city in cities:
        xy = m(city[1], city[0])
        x = xy[0]
        y = xy[1]
        m.scatter(x, y, marker='o', s=markersize, edgecolor='k', c='w',
                  linewidth=2, zorder=zorder)
        if names:
            plt.text(x + 0.02 * (m.xmax - m.xmin), y + 0.02 * (m.ymax - m.ymin), 
                     s = city[2], fontsize=fontsize, zorder=zorder,
                     fontweight='roman',
                     bbox=dict(facecolor='w', fill=True, edgecolor='w', 
                               alpha=0.75, boxstyle='round')
                     )


def plot_beachballs(m, cat, **kwargs):
    """
    plot moment tensors as beachballs, color by depth if requested
    """
    depths, mags = [], []
    for event in cat:
        depths.append(event.preferred_origin().depth * 1E-3)
        mags.append(event.preferred_magnitude().mag)
    
    # normalize the colormap, create a discrete colorbar
    colormap = plt.cm.plasma
    max_depth = myround(max(depths), 5, 'up')
    normalize = mpl.colors.BoundaryNorm(range(0, max_depth, 5), colormap.N)

    # normalize the magnitudes
    magnitudes = normalize_a_to_b(mags, a=1.2E4, b=2.1E4) 

    for i, event in enumerate(cat):
        event_beachball(m, event, fm_type="strike_dip_rake", 
                        facecolor=colormap(normalize(depths[i])), 
                        width=magnitudes[i],
                        **kwargs
                        )
    
    # create a colormap for the colorbar
    sm = plt.cm.ScalarMappable(cmap=colormap, norm=normalize)
    sm.set_array([])
    cbar = plt.colorbar(sm, extend='max', shrink=0.35, aspect=10)
    cbar.set_label("depth (km)", rotation=270, labelpad=15, fontsize=15)
    cbar.ax.invert_yaxis()
    cbar.ax.tick_params(axis='y', direction='in', length=0, width=1.75)


def plot_interface(m, fid, **kwargs):
    """
    plot the plate interface 
    """
    zorder = kwargs.get("zorder", 100)
    colormap = kwargs.get("colormap", None)

    plate = np.load(fid)
    lats = plate[:, 0]
    lons = plate[:, 1]
    depths = plate[:, 2]

    x, y = m(lons, lats)
    xi, yi = np.mgrid[min(x):max(y):1000,
                      min(x):max(y):2000]
    zi = griddata((x, y), depths, (xi, yi), 'cubic')

    import ipdb;ipdb.set_trace()

    cs = m.contour(xi, yi, zi, vmin=0, zorder=zorder, cmap=plt.cm.plasma)

    # m.scatter(x, y, c=depths, cmap=colormap)


def srcrcv_map():
    """
    main
    """
    # Plotting kwargs
    figsize = (10, 12)
    dpi = 100
    coastline_linewidth = 2.25
    continent_color= "w"
    lake_color = "w"
    fill_color = "w"
    station_markersize = 110
    station_zorder = 100
    beachball_linewidths = 1.6
    beachball_zorder = station_zorder - 1
    city_zorder = station_zorder + 1
    city_markersize = 100
    city_fontsize = 12
    city_names = False
    scalebar_location = "lower-right"
    plate_zorder = station_zorder - 2

    # User Parameters
    # Lat_min: 20 km south of KHZ (173.359, -42.416)
    # Lat_max: 20 km north of KAZ (175.1612, -37.1041)
    # Lon_min: 20 km west of NNZ (173.3795, -41.2171)
    # Lon_max: 20 km east of MXZ (179.3066, -37.5623)
    map_corners={'lat_min': -42.6, 'lat_max': -36.925,
                 'lon_min': 173.14, 'lon_max': 178.5}

    # File paths
    cat = read_events("./charlie_trial.xml")
    stations_fid = "./TRIALS_STATIONS_78"
    plate_fid = "./plate_interface.npy"


    # Workflow
    f = plt.figure(figsize=figsize, dpi=dpi)
    ax = f.add_subplot(111)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(2)

    m = initiate_basemap(map_corners, continent_color=continent_color, 
                         lake_color=lake_color, fill_color=fill_color,
                         coastline_linewidth=coastline_linewidth, scalebar=True,
                         scalebar_location=scalebar_location
                         )

    plot_stations(m, stations_fid, markersize=station_markersize,
                  zorder=station_zorder)
    plot_beachballs(m, cat, zorder=beachball_zorder,
                    linewidth=beachball_linewidths)
    plot_cities(m, zorder=city_zorder, markersize=city_markersize,
                fontsize=city_fontsize, names=city_names)
    plot_interface(m, plate_fid, zorder=plate_zorder)
    # plot_geonet_active_faults(m)

    plt.show()
    

if __name__ == "__main__":
    srcrcv_map()
