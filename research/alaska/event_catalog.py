"""
The catalog used for initial testing purposes for Pyatoa + Seisflows
26.7.19
"""
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime, Catalog, read_events
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn.header import FDSNNoDataException


def cut_cat(cat, indices=None, method="remove"):
    """
    Obspy's catalog object has no delete, so here it is!

    :type cat: obspy.Catalog
    :param cat: catalog object to cut down
    :type indices: list of ints
    :param indices: indices of the catalog to remove
    :type method: str
    :param method: "remove" to delete 'indices' from 'cat', "keep" to retain
    :rtype: obspy.Catalog
    :return: new catalog object cut from old one
    """
    if indices is None:
        return cat

    new_cat = Catalog()
    if method == "remove":
        events = [e for i, e in enumerate(cat) if i not in indices]
    elif method == "keep":
        events = [e for i, e in enumerate(cat) if i in indices]
    new_cat.extend(events)

    return new_cat


def plot_cat(cat, evnts, tag=""):
    """
    Plot the catalog with a custom title based on the

    :type cat: obspy.Catalog
    :param cat: catalog to plot
    :type evnts: dictionary
    :param evnts: parameter dictionary from function event_trials()
    :type tag: str
    :param tag: tag for saving the plot name
    """
    # Plot the catalog
    evnts['starttime'].precision = 1
    evnts['endtime'].precision = 1

    title = "\n".join([
        "{c} {n} events".format(c=evnts['name'], n=len(cat)),
        "{s} - {e}".format(s=evnts['starttime'], e=evnts['endtime']),
        "M{a} - {b}".format(a=evnts['minmagnitude'], b=evnts['maxmagnitude']),
        "LLC: {latmin},{lonmin} / URC: {latmax},{lonmax}\n".format(
            latmin=evnts['minlatitude'], lonmin=evnts['minlongitude'],
            latmax=evnts['maxlatitude'], lonmax=evnts['maxlongitude'])
    ])

    cat.plot(projection='local', resolution='l', continent_fill_color='w',
             water_fill_color='w', color='date',
             outfile="./{}{}".format(evnts['name'], tag), title=title
             )

    plt.close('all')


def event_trials(choice, sep_km=None, desired_length=None):
    """
    For a given set of trial parameters, retrieve an event catalog and then
    slim it down with various check parameters; make sure events are not too 
    closely grouped, and that they contain moment tensor information from GeoNet
    Plot the final catalog, make beachballs for plotting purposes, return
    an .xml catalog

    :type choice: str
    :param choice: choice of parameter trial
    :type csv_file: str
    :param csv_file: path to GeoNet csv file with moment tensor information
    :type desired_length: int
    :param desired_length: desired length of catalog, event removal will try to 
    :return:
    """
    trials = {
        "nalaska_eq_Mgt2": {
            "name": "nalaska_earthquakes_Mgt2",
            "starttime": UTCDateTime("1990-01-01T00:00:00"), 
            "endtime": UTCDateTime("2023-01-01T00:00:00"),
            "minmagnitude": 2.,
            "maxmagnitude": 9.,
            "minlatitude": 64,  # LLC
            "minlongitude": -169,  # LLC
            "maxlatitude": 72,  # URC
            "maxlongitude": -137,  # URC
            "maxdepth": 50.,
            "client": "IRIS",
          }
    }

    evnts = trials[choice]
    
    # Create the catalog and plot the raw event catalog
    original_cat_fid = f"{cat_name}_original.xml"
    if not os.path.exists(original_cat_fid):
        c = Client(evnts["client"])
        original_cat = c.get_events(starttime=evnts['starttime'],
                                    endtime=evnts['endtime'],
                                    minmagnitude=evnts['minmagnitude'],
                                    maxmagnitude=evnts['maxmagnitude'],
                                    minlatitude=evnts['minlatitude'],
                                    minlongitude=evnts['minlongitude'],
                                    maxlatitude=evnts['maxlatitude'],
                                    maxlongitude=evnts['maxlongitude'],
                                    maxdepth=evnts['maxdepth'],
                                    orderby="magnitude-asc"
                                    )
        print("catalog has {} events".format(len(original_cat)))
        # Delete unncessary attributes to save space
        for attr in ["picks", "amplitudes", "station_magnitudes"]:
            for event in original_cat:
                try:
                    delattr(event, attr)
                except KeyError:
                    continue
        original_cat.write(original_cat_fid, format="QUAKEML")


if __name__ == "__main__":
    cat_name = "nalaska_eq_Mgt2"
    event_trials(cat_name)

