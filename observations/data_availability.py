"""
A good way to show what days have data for a given network based on
availability of mseed files
"""
import os
import glob
import numpy as np

path = "./"
years = glob.glob(os.path.join(path, "????"))
years.sort()
for year in years:
    netstas = []
    filled_years = np.array()
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
        netsas.append(netstaname)
        filled_years = np.append(filled_years, station_full_year)
        import ipdb;ipdb.set_trace()




