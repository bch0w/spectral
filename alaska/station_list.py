"""
Pyatoa and SPECFEM require station lists for various reasons. The most important
information for a station are:
    code, network, latitude, longitude, response
pair that describes location. Elevation is not necessarily important, but is
useful to have. Pyatoa requires station start and endtimes to determine
data availability as well as for map plotting.

This script creates a master station list for the north island New Zealand
adjoint tomography problem. Each function gathers relevant station information
and assembles it into an obspy.network class; these are then combined into
an obspy.inventory class and written out as a StationXML file. Functionality
also allows for writing out to SPECFEM text file format if necessary.
"""
import os
import numpy as np
from obspy.clients.fdsn import Client
from obspy.clients.nrl import NRL
from obspy import Inventory, UTCDateTime
from obspy.core.inventory.network import Network
from obspy.core.inventory.station import Station
from obspy.core.inventory.channel import Channel
from obspy.core.inventory.util import Site

# from pyatoa.utils.srcrcv import merge_inventories


def get_from_iris(level="station"):
    """
    Return available broadband stations from IRIS

    :type level: str
    :param level: level to propogate network creation
    :rtype: obspy.core.inventory.network.Network
    """
    c = Client("IRIS")

    # A whole-Alaska bounding box including western Canada and Aleutians
    lat_min = 64.5  # South
    lat_max = 72.  # North
    lon_min = -168.    # West
    lon_max = -140.  # East

    starttime = UTCDateTime("2000-01-01")
    networks = ["AK",  # Alaska Regional Network
                "AT",  # National Tsunami Warning Center (Alaska Seismic Network)
                "AV",  # Alaska Volcano Observatory
                "CN",  # Canadian National Seismograph Network
                "II",  # Global Seismograph Network (GSN)
                "IU",  # Global Seismograph Network (GSN)
                "NY",  # Yukon-Northwest Seismic Network
                "TA",  # USArray Transportable Array
                "US",  # United States National Seismic Network
                "YO",  # Yukon Observatory
                ]
    channels = ["BH?", "BL?", "HH?", "HL?"]

    # ObsPy doesnt accept lists or bracketed wildcards, must be comma separated
    networks = ",".join(networks)
    channels = ",".join(channels)

    inv_alaska = c.get_stations(network=networks, station="*", channel=channels,
                                minlatitude=lat_min, maxlatitude=lat_max,
                                minlongitude=lon_min, maxlongitude=lon_max,
                                level=level, starttime=starttime,
                                endtime=UTCDateTime()
                                )

    return inv_alaska


def export_specfem(inv, filename="master_station_list.txt"):
    """
    Specfem requires station information in the form

    STATION NETWORK LATITUDE LONGITUDE ELEVATION BURIAL

    We can leave elevation and burial as 0.0 because our topography data takes
    care of elevation (stations assumed to be on the surface, which is
    approximately true. Take a stationxml and return a text file (no formatter)
    that contains the correct information

    :return:
    """
    template = ("{station:>6}{network:>6}    {latitude:6.4f}    "
                "{longitude:6.4f}    0.0    0.0\n"
                )
    with open(filename, 'w') as f:
        for net in inv:
            for sta in net:
                f.write(template.format(station=sta.code, network=net.code,
                                        latitude=sta.latitude,
                                        longitude=sta.longitude)
                        )


def export_seed_fmt(inv):
    """
    Pyatoa requires a "seed" convention directory structure containing
    individual RESPONSE xml files. Create that here based on Obspy inventory

    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory object to export
    :return:
    """
    path = "./seed/RESPONSE"
    dir_structure = '{sta}.{net}'
    file_template = 'RESP.{net}.{sta}.{loc}.{cha}'

    if not os.path.exists(path):
        os.makedirs(path)

    # Parse through inventory
    for net in inv:
        for sta in net:
            # Create the directory to store individual channels
            sta_dir = os.path.join(path, dir_structure.format(
                sta=sta.code, net=net.code)
                                   )
            if not os.path.exists(sta_dir):
                os.makedirs(sta_dir)
            for cha in sta:
                # Select the channel as a
                channel_as_inv = inv.select(network=net.code, station=sta.code,
                                            channel=cha.code
                                            )
                # If select returns properly
                if channel_as_inv:
                    fid_out = os.path.join(
                        path, dir_structure.format(sta=sta.code, net=net.code),
                        file_template.format(net=net.code, sta=sta.code,
                                             loc=cha.location_code,
                                             cha=cha.code)
                    )
                    channel_as_inv.write(fid_out, format='STATIONXML')
                else:
                    print("{}.{} could not be selected".format(sta.code,
                                                               cha.code)
                          )


if __name__ == "__main__":
    # Parameters
    level = "channel"  # channel, station
    write_to = "alaska_bb_stations.xml"  # Name of output .xml file
    export_to_specfem = True
    export_to_seed_fmt = False
    plot = True

    # Create the Inventory
    master_inventory = get_from_iris(level=level)


    # Export to various output formats
    if write_to:
        master_inventory.write(write_to, format="STATIONXML")
    if export_to_specfem:
        export_specfem(master_inventory)
    if export_to_seed_fmt:
        export_seed_fmt(master_inventory)
    if plot:
        master_inventory.plot(projection="local", resolution="l",
                              continent_fill_color="w", water_fill_color="w",
                              label=True, color_per_network=True,
                              outfile="./nalaska.png"
                              )
