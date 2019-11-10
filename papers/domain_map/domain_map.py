"""
Overly-detailed publication quality source receiver map
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from obspy import read_inventory, read_events
from pyatoa.utils.visuals.mapping import standalone_map


def plot_inv(m, inv):
    """
    plot stations from an obspy inventory object
    """
    # place changeable variables up high for quick change
    markersize = 100
    zorder = 100

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

def plot_stations(m, stations_fid):
    """
    plot stations from a Specfem STATION file
    """
    # place changeable variables up high for quick change
    markersize = 100
    zorder = 100

    station_info = np.loadtxt(stations_fid, usecols=(0,1,2,3), dtype=str)
    x, y = [], []
    for station in station_info:
        sta, net, lat, lon = station
        name = f"{net}.{sta}"
        xy = m(float(lon), float(lat))
        x.append(xy[0])
        y.append(xy[1])

    m.scatter(x, y, marker=11, s=markersize, edgecolor='k', c='k',
              linestyle='-', linewidth=1.25, cmap='viridis', 
              zorder=zorder)


def domain_map():
    """
    main
    """
    # Parameters
    # inv = read_inventory("./master_inventory.xml")
    cat = read_events("./charlie_trial.xml")
    stations_fid = "./TRIALS_STATIONS_78"

    map_corners={'lat_min': -42.5007, 'lat_max': -36.9488,
                 'lon_min': 172.9998, 'lon_max': 179.5077}

    # workflow
    f, m = standalone_map(map_corners, show="hold")
    plot_stations(m, stations_fid)

    plt.show()

    

if __name__ == "__main__":
    domain_map()
