"""
A good way to show what days have data for a given network based on
availability of mseed files
"""
import os
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

mpl.rcParams['font.size'] = 12

def myround(x, base=5, choice='near'):
    """
    Round value x to nearest base, round 'up','down' or to 'near'est base

    :type x: float
    :param x: value to be rounded
    :type base: int
    :param base: nearest integer to be rounded to
    :type choice: str
    :param choice: method of rounding, 'up', 'down' or 'near'
    :rtype roundout: int
    :return: rounded value
    """
    if choice == 'near':
        roundout = int(base * round(float(x)/base))
    elif choice == 'down':
        roundout = int(base * np.floor(float(x)/base))
    elif choice == 'up':
        roundout = int(base * np.ceil(float(x)/base))

    return roundout


# define service dates


path = "./"
years = glob.glob(os.path.join(path, "????"))
years.sort()
for year in years:
    earliest, latest = 1E3, 0
    netstas = []
    filled_years = None
    networks = glob.glob(os.path.join(year, "*"))
    networks.sort()
    for network in networks:
        stations = glob.glob(os.path.join(network, "*"))
        stations.sort()
        for station in stations:
            channels = glob.glob(os.path.join(station, '*'))
            channels.sort()
            station_full_year = np.zeros(366)
            for channel in channels:
                fileids = glob.glob(os.path.join(channel, "*"))
                fileids.sort()
                # File ids will add to the bit in the full year
                for fileid in fileids:
                    jday = int(fileid.split('.')[-1])
                    # determine what the start and end dates are
                    if jday < earliest:
                        earliest = jday
                    elif jday > latest:
                        latest = jday
                    station_full_year[jday] += 1

            netstaname = "{}".format(os.path.basename(station))
            netstas.append(netstaname)
            if filled_years is None:
                filled_years = station_full_year
            else:
                filled_years = np.vstack([filled_years, station_full_year])

        # reduce data to the available bound
        filled_years = filled_years[:, earliest:latest]

        fig, ax = plt.subplots(figsize=(20,10)) 
        im = ax.imshow(filled_years, extent=(earliest, latest, len(netstas), 0),
                       aspect="auto", cmap="Oranges")
        year_ = os.path.basename(year)
        ax.set_title(f"{year_} BEACON DATA AVAILABILITY")
        ax.set_xlabel(f"Julian Day ({year_})")
        ax.set_ylabel("Station Code")

        # colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='2.5%', pad=0.05)
        cbar = fig.colorbar(im, cax=cax, orientation='vertical', 
                            values=[0,1,2,3], ticks=[0,1,2,3])
        cbar.set_label('number of available channels', rotation=90, 
                       labelpad=-15)
        
        # grid
        ax.grid(True, which='major', linestyle='-', c="k", linewidth=1.25)
        ax.grid(True, which='minor', linestyle='-', c="k", linewidth=0.5, 
                alpha=0.5)
        ax.set_xticks(np.arange(myround(earliest, 5, 'near'), latest, 5))
        ax.set_xticks(np.arange(myround(earliest, 5, 'near'), latest, 1), 
                      minor=True)
        ax.set_yticks(np.arange(len(netstas)))
        ax.set_yticklabels(netstas)
        plt.setp(ax.get_xticklabels(), rotation=45, fontsize=13)
        plt.setp(ax.get_yticklabels(), rotation=45, fontsize=13)

        plt.show()
        # plt.savefig("data_{}.png".format(os.path.basename(year)))






