"""
A good way to show what days have data for a given network based on
availability of mseed files
"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt


path = "./"
years = glob.glob(os.path.join(path, "????"))
years.sort()
for year in years:
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
                    station_full_year[jday] += 1

            netstaname = "{}.{}".format(os.path.basename(network),
                                        os.path.basename(station))
            netstas.append(netstaname)
            if filled_years is None:
                filled_years = station_full_year
            else:
                filled_years = np.vstack([filled_years, station_full_year])

        fig, ax = plt.subplots(figsize=(10,10)) 
        im = ax.imshow(filled_years, aspect="auto")
        ax.set_xticks(np.arange(len(station_full_year)))
        ax.set_yticks(np.arange(len(netstas)))
        ax.set_yticklabels(netstas)
        plt.setp(ax.get_xticklabels(), rotation=90)
        plt.grid()
        ax.set_title("{} DATA".format(os.path.basename(year)))
        plt.savefig("data_{}.png".format(os.path.basename(year)))






