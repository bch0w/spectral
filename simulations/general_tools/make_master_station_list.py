"""
Pyatoa and SPECFEM require station lists for various reasons. The most important
information for a station is the code, related network and latitude longitude
pair that describes location. Elevation is not necessarily important, but is
useful to have. Pyatoa requires station start and endtimes to determine
data availability as well as for map plotting.

This script creates a master station list for the north island New Zealand
adjoint tomography problem. Each function gathers relevant station information
and assembles it into an obspy.network class; these are then combined into
an obspy.inventory class and written out as a StationXML file. Functionality
also allows for writing out to SPECFEM text file format if necessary.
"""

import glob

import numpy as np
from obspy.clients.fdsn import Client
from obspy import Inventory, UTCDateTime, read_inventory
from obspy.core.inventory.network import Network
from obspy.core.inventory.station import Station
from obspy.core.inventory.util import Site


def fetch_geonet_network():
    """
    return geonet broadband stations in an obspy network, no response
    TO DO: manual ignore some TVZ stations to keep clutter down in that area

    :return:
    """
    c = Client("GEONET")
    lat_lon = [-42.8, -37.1, 173, 179.4]  # srtm30p 550 641 mesh

    # GeoNet Broadband instruments
    inv_geonet = c.get_stations(network='NZ', station='*Z', channel='HH?',
                                minlatitude=lat_lon[0], maxlatitude=lat_lon[1],
                                minlongitude=lat_lon[2],
                                maxlongitude=lat_lon[3],
                                level="station",
                                starttime=UTCDateTime('2000-01-01'),
                                endtime=UTCDateTime()
                                )

    return inv_geonet[0]


def fetch_hobitss_network():
    """
    return hobitss broadband OBS sensors (trillium compacts)
    :return: 
    """
    c = Client("IRIS")
    inv_hobitss = c.get_stations(network='YH', station="LOBS*", location='',
                                 level="station"
                                 )
    return inv_hobitss[0]


def build_sahke_network():
    """
    sahke transect broadband stations aren't fetchable through fdsn webservies,
    build a network object from the available information
    :return:
    """

    # station: start, stop, lat, lon
    network_code = "X2"
    station_info = [["LE4", "2010-136", "2010-331", -41.3579, 175.6919],
                    ["LTN6", "2010-193", "2010-349", -41.1033, 175.3238],
                    ["T004", "2010-088", "2010-255", -41.3403, 175.6688],
                    ["T007", "2010-041", "2010-123", -41.3041, 175.6513],
                    ["T010", "2010-135", "2010-348", -41.2520, 175.5825],
                    ["T014", "2010-034", "2010-350", -41.2075, 175.5063],
                    ["T016", "2010-088", "2010-322", -41.1893, 175.4737],
                    ["T018", "2010-055", "2010-349", -41.1715, 175.3850],
                    ["T020", "2010-089", "2010-261", -41.1251, 175.3497]
                    ]
    default_elevation = 123456.0
    default_site = Site(name="SAHKE")

    stations = []
    for stalist in station_info:
        station_ = Station(code=stalist[0], start_date=UTCDateTime(stalist[1]),
                           end_date=UTCDateTime(stalist[2]),
                           latitude=stalist[3], longitude=stalist[4],
                           elevation=default_elevation,
                           site=default_site, creation_date=UTCDateTime()
                           )
        stations.append(station_)
    network = Network(code=network_code,
                      description="SAHKE2 Broadband Transect",
                      stations=stations
                      )
    return network


def build_bannister_network():
    """
    bannister broadband data not available through fdsn
    :return:
    """
    network_code = "X1"
    # station, lat, lon, depth, start, end
    station_info = [['GA01', -39.0331, 177.8549, 33.000, 2011305, 9999001],
                    ['GA02', -38.7793, 177.8695, 220.000, 2011305, 9999001],
                    ['GA03', -38.6059, 177.9959, 46.000, 2011305, 9999001],
                    ['GA04', -38.5529, 178.1548, 210.000, 2011305, 9999001],
                    ['GA05', -38.5844, 177.8048, 85.000, 2011305, 9999001],
                    ['GA06', -39.1110, 177.9161, 159.000, 2012001, 9999001],
                    ['GA07', -38.4322, 178.0478, 569.000, 2012001, 9999001],
                    ['GA08', -38.8312, 177.7095, 358.000, 2012001, 9999001],
                    ['GA09', -38.5194, 177.9363, 107.000, 2013152, 9999001],
                    ['GA10', -38.9087, 177.4627, 57.000, 2013152, 9999001],
                    ['GA11', -38.4932, 177.6711, 224.000, 2014060, 9999001],
                    ['HD01', -38.4802, 175.9473, 470.000, 2009244, 2010091],
                    ['HD02', -38.6275, 175.9196, 600.000, 2009244, 2012306],
                    ['HD03', -38.5497, 176.0564, 630.000, 2009244, 2010091],
                    ['HD04', -38.4930, 176.2204, 440.000, 2009244, 9999001],
                    ['HD05', -38.4663, 176.2627, 445.000, 2009244, 2011091],
                    ['HD06', -38.3932, 176.0619, 280.000, 2009244, 2011091],
                    ['HD07', -38.3214, 176.1626, 480.000, 2009244, 2010091],
                    ['HD08', -38.6300, 176.3063, 520.000, 2009244, 2011091],
                    ['HD09', -38.6675, 176.1798, 450.000, 2010001, 2010305],
                    ['HD10', -38.5482, 176.3669, 390.000, 2010001, 2011091],
                    ['HD11', -38.6320, 176.2606, 320.000, 2010091, 2011305],
                    ['HD12', -38.3711, 176.1570, 360.000, 2010091, 2011091],
                    ['HD13', -38.4546, 176.3458, 340.000, 2010091, 2011091],
                    ['HD14', -38.4594, 176.1714, 340.000, 2010091, 2010305],
                    ['HD15', -38.4872, 176.0043, 660.000, 2010091, 2011091],
                    ['HD16', -38.4408, 175.9444, 420.000, 2010091, 2011091],
                    ['HD17', -38.5483, 175.0556, 600.000, 2010091, 2011091],
                    ['HD18', -38.5283, 176.4573, 375.000, 2010305, 2011091],
                    ['HD19', -38.3762, 176.3696, 325.000, 2010305, 2011091],
                    ['HD20', -38.5842, 176.1467, 430.000, 2010305, 2011091],
                    ['HD21', -38.5052, 176.0877, 515.000, 2010305, 2011091],
                    ['HD22', -38.5666, 176.1907, 385.000, 2010305, 2011001],
                    ['HD23', -38.6253, 175.9533, 520.000, 2010305, 2011091],
                    ['HD24', -38.5671, 176.0895, 560.000, 2010305, 2011091],
                    ['HD25', -38.5855, 176.2945, 305.000, 2010305, 2011091],
                    ['HD26', -38.5696, 175.9547, 610.000, 2010305, 2011091],
                    ['HD27', -38.4324, 176.2374, 530.000, 2010305, 2011091],
                    ['HD28', -38.4424, 176.1580, 380.000, 2010305, 2011091],
                    ['HD29', -38.4012, 176.1991, 460.000, 2010305, 2011091],
                    ['HD30', -38.5291, 175.9377, 480.000, 2010305, 2011091],
                    ['HD31', -38.3019, 176.3053, 525.000, 2010305, 2011091],
                    ['HD32', -38.7035, 176.1434, 558.000, 2010335, 2011091],
                    ['HD33', -38.5031, 176.2634, 305.000, 2011001, 2011091],
                    ['HD34', -38.4764, 176.0497, 515.000, 2011032, 2011091],
                    ['HD35', -38.5331, 176.0008, 500.000, 2011032, 2011060],
                    ['HD36', -38.4912, 176.1677, 340.000, 2011032, 2011091],
                    ['HD37', -38.4177, 176.1366, 355.000, 2011032, 2011091],
                    ['HD38', -38.5190, 176.2041, 325.000, 2011032, 2011091],
                    ['HD39', -38.4868, 176.3221, 295.000, 2011032, 2011091],
                    ['HD50', -38.3396, 176.2687, 361.000, 2015305, 9999001],
                    ['HD51', -38.2435, 176.2483, 380.000, 2015305, 9999001],
                    ['HD53', -38.2659, 176.4812, 362.000, 2015335, 9999001],
                    ['HD54', -38.2815, 176.3781, 480.000, 2015335, 9999001],
                    ['HD55', -38.2815, 176.5605, 469.000, 2015335, 9999001],
                    ['HD56', -38.3184, 176.5538, 439.000, 2015335, 9999001],
                    ['HD57', -38.1872, 176.3714, 520.000, 2015335, 9999001],
                    ['HD58', -38.3019, 176.3053, 492.000, 2015335, 9999001],
                    ['HD59', -38.2815, 176.3781, 440.000, 2015335, 9999001],
                    ['HD60', -38.2863, 176.1992, 389.000, 2016061, 9999001]
                    ]

    stations = []
    default_site = Site(name="Banni")
    for stalist in station_info:
        if stalist[5] == 9999001:
            end_date = UTCDateTime()
        else:
            end_date = UTCDateTime(stalist[5])
        station_ = Station(code=stalist[0], latitude=stalist[1],
                           longitude=stalist[2], elevation=stalist[3],
                           start_date=UTCDateTime(str(stalist[4])),
                           end_date=end_date, site=default_site,
                           creation_date=UTCDateTime()
                           )
        stations.append(station_)
    network = Network(code=network_code,
                      description="Stephen Bannister's Temporary Network",
                      stations=stations
                      )
    return network


def read_fathom_net():
    """
    get fathom network from dataless file, don't include channel information
    :return:
    """
    network_code = "XX"
    path_ = ("/Users/chowbr/Documents/subduction/mseeds/FATHOM/DATALESS/"
             "XX.RDF.DATALESS")
    inv = read_inventory(path_)

    stations = []
    for sta_ in inv[0]:
        station_ = Station(code=sta_.code, latitude=sta_.latitude,
                           longitude=sta_.longitude, elevation=sta_.longitude,
                           start_date=sta_.start_date, end_date=sta_.end_date,
                           site=sta_.site, creation_date=UTCDateTime()
                           )
        stations.append(station_)
    network = Network(code=network_code,
                      description="FATHOM Network",
                      stations=stations
                      )
    return network


def generate_master_list():
    """
    call all functions in this script to create a master list of station data
    :return:
    """
    geonet = fetch_geonet_network()
    hobitss = fetch_hobitss_network()
    sahke = build_sahke_network()
    bannister = build_bannister_network()
    fathom = read_fathom_net()

    return Inventory(networks=[geonet, hobitss, sahke, bannister, fathom],
                     source="homemade"
                     )


def output_inventory_to_specfem_format(inventory):
    """
    specfem requires station information in the form
    STATION NETWORK LATITUDE LONGITUDE ELEVATION BURIAL
    we can leave elevation and burial as 0.0 because our topography data takes
    care of elevation (stations assumed to be on the surface, which is
    approximately true. Take a stationxml and return a text file (no formatter)
    that contains the correct information
    :return:
    """
    template = ("{station:>6}{network:>6}    {latitude:.4f}"
                "{longitude:.4f}    0.0    0.0\n"
                )
    with open('master_station_list.txt','w') as f:
        for net in inv:
            for sta in net:
                f.write(template.format(station=sta.code, network=net.code,
                                        latitude=sta.latitude,
                                        longitude=sta.longitude)
                        )


def output_inventory_to_numpy_array(inventory):
    """
    .npz files are much more compact and easier to read/write than stationxml's
    which contain a lot of fluff. this function will parse an inventory object
    and return a .npz file with the necessary information for pyatoa
    :param inventory:
    :return:
    """
    # to do, lol


if __name__ == "__main__":
    master_inventory = generate_master_list()
    master_inventory.write('master_inventory.xml', format='STATIONXML')
