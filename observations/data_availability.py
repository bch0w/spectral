"""
A good way to show what days have data for a given network based on
availability of mseed files
/"""
import os
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime, read_events
from mpl_toolkits.axes_grid1 import make_axes_locatable


mpl.rcParams['font.size'] = 15


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

def plot_availability(path="./", cat=None):
    """
    Plot data availability in a bar chart by year for easy visual identification

    :type path: str
    :param path: path to data in seed structure
    :type cat: obspy.Catalog
    :param cat: If a catalog object is given, event dates will be highlighted
        for easy event identification
    """
    years = glob.glob(os.path.join(path, "????"))
    years.sort()
    for year in years:
        earliest, latest = 1E3, 0
        netstas = []
        filled_years = None
        # networks = glob.glob(os.path.join(year, "*"))
        networks = glob.glob(os.path.join(year, "XX"))
        if not networks:        
            continue
        networks.sort()
        for network in networks:
            # stations = glob.glob(os.path.join(network, "*"))
            # stations.sort()
            stations = [os.path.join(network, f"RD{_:0>2}") for _ in range(1, 23)]
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

                netstaname = f"{os.path.basename(station)}"
                netstas.append(netstaname)
                if filled_years is None:
                    filled_years = station_full_year
                else:
                    filled_years = np.vstack([filled_years, station_full_year])

        np.savez(f"{year.split('/')[1]}", filled_years) 
        if "2019" in year:
            a = 1/0 
        else:
            continue

        # reduce data to the available bound
        # filled_years = filled_years[:, earliest:latest]
        cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["white", 
                                                                 "red",
                                                                 "yellow",
                                                                 "green"])
        fig, ax = plt.subplots(figsize=(12,4),dpi=400) 
        im = ax.imshow(filled_years, 
                       extent=(1, 365, len(netstas), 0),
                       # extent=(earliest, latest, len(netstas), 0),
                       aspect="auto", cmap=cmap) #aspect='auto'

        year_ = os.path.basename(year)
        ax.set_title(f"{year_}")
        ax.set_xlabel(f"Julian Day ({year_})")
        ax.set_ylabel("Station Code")

        # colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='2.5%', pad=0.05)
        cbar = fig.colorbar(im, cax=cax, orientation='vertical', 
                            values=[0,1,2,3], ticks=[0,1,2,3])
        cbar.set_label('number of available channels', rotation=90, 
                       labelpad=-12, fontsize=12)
        
        # grid
        ax.grid(True, which='major', linestyle='-', c="k", linewidth=.6, 
                alpha=.5)
        ax.grid(True, which='minor', linestyle='-', c="k", linewidth=0.3, 
                alpha=0.5)
        ax.set_xticks(np.arange(-.5, 365.5, 10))
        ax.set_xticks(np.arange(-.5, 365.5, 1), minor=True)
        ax.set_xticklabels(labels=[1] + [str(_) for _ in np.arange(10, 366, 10)])
        # ax.set_xticks(np.arange(myround(earliest, 10, 'near'), latest, 10))
        # ax.set_xticks(np.arange(myround(earliest, 10, 'near'), latest, 1), 
        #               minor=True)
        ax.set_yticks(np.arange(len(netstas)))
        ax.set_yticklabels(netstas)
        plt.setp(ax.get_xticklabels(), rotation=45, fontsize=13, ha="center")
        plt.setp(ax.get_yticklabels(), rotation=0, fontsize=13, va="top")

        
        # plot catalog
        if cat:
            plot_catalog(ax, cat, int(os.path.basename(year)))
        # plt.show()
        plt.savefig(f"data_{os.path.basename(year)}.png")


def plot_catalog(ax, cat, year):
    """
    Plot event dates in a catalog as vertical lines
        
    :type ax: matplotlib.axis
    :param ax: axis to plot to
    :type cat: obspy.Catalog
    :param cat: catalog to check
    :type year: int 
    :param year: year to check
    """
    days = []
    for event in cat:
        event_time = event.preferred_origin().time

        # Make sure the event falls into the time frame of the data
        startday, endday = ax.get_xbound()
        checkday = (event_time.julday > startday and event_time.julday < endday)
        
        if event_time.year == year and checkday:
            days.append(event_time.julday)

            # Set the y bound incase multiple events on the same day
            ymin, ymax = ax.get_ybound()
            y_text = (ymax - ymin) / (days.count(event_time.julday) + 1)

            ax.axvline(x=event_time.julday, ymin=0, ymax=1, color="k",
                       linestyle="--", linewidth=2)
            ax.text(x=event_time.julday, 
                    y=(ax.get_ybound()[1] - ax.get_ybound()[0]) / 2,
                    s=event.resource_id.id.split('/')[1], 
                    color="k", rotation=90)


def check_availability(origintime, path="./",  return_time=False):
    """
    Check which stations have data for a given day, can also check how much data
    is available for a given day

    :type origintime: str
    :param origintime: datetime string of event origin, e.g. 2000-01-01T00:00:00
    :type path: str
    :param path: path to data in seed format
    :type return_time: bool
    :param return_time: also print the timing availability, takes longer as this
        requires reading in the mseeds using obspy
    """
    print(origintime)
    origin = UTCDateTime(origintime)

    years = glob.glob(os.path.join(path, "????"))
    years.sort()
    for year in years:
        if int(os.path.basename(year)) != origin.year:
            continue
        networks = glob.glob(os.path.join(year, "*"))
        networks.sort()
        for network in networks:
            stations = glob.glob(os.path.join(network, "*"))
            stations.sort()
            for station in stations:
                channels = glob.glob(os.path.join(station, '*'))
                channels.sort()
                for channel in channels:
                    fileids = glob.glob(os.path.join(channel, "*"))
                    fileids.sort()
                    for fileid in fileids:
                        jday = int(fileid.split('.')[-1])
                        if jday == origin.julday:
                            if return_time:
                                st_ = read(fileid)
                                print(f"\t{os.path.basename(fileid)} "
                                      f"{st_[0].stats.starttime} "
                                      f"{st_[0].stats.endtime}")
                            else:
                                print(f"\t{os.path.basename(fileid)}")
        

if __name__ == "__main__":
    # cat = read_events("./fullscale_w_mt.xml")
    cat = None
    plot_availability("../", cat=cat)

