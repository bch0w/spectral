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
from obspy.core.inventory.util import Site, Equipment

from pyatoa.utils.srcrcv import merge_inventories


def geonet_north_island(level="station"):
    """
    Return geonet broadband stations in an obspy network

    :type level: str
    :param level: level to propogate network creation
    :rtype: obspy.core.inventory.network.Network
    """
    c = Client("GEONET")

    # Extended North Island coverage, Kaikoura to East Cape
    lat_min = -42.8
    lat_max = -37.1
    lon_min = 173
    lon_max = 179.4

    starttime = UTCDateTime("2000-01-01")

    # GeoNet Broadband instruments with standard station code Z
    inv_geonet = c.get_stations(network="NZ", station="*Z", channel="HH?",
                                minlatitude=lat_min, maxlatitude=lat_max,
                                minlongitude=lon_min, maxlongitude=lon_max,
                                level=level, starttime=starttime,
                                endtime=UTCDateTime()
                                )

    inv_hold = c.get_stations(network="NZ", station="*Z", channel="BH?",
                              minlatitude=lat_min, maxlatitude=lat_max,
                              minlongitude=lon_min, maxlongitude=lon_max,
                              level=level, starttime=starttime,
                              endtime=UTCDateTime()
                              )
    inv_geonet = merge_inventories(inv_geonet, inv_hold)
    # Stations with unique names (Wellington and Baring Head)
    for station in ["WEL", "BHW"]:
        inv_hold = c.get_stations(network="NZ", station=station, 
                                  channel="HH?", minlatitude=lat_min, 
                                  maxlatitude=lat_max, minlongitude=lon_min, 
                                  maxlongitude=lon_max, level=level, 
                                  starttime=starttime, endtime=UTCDateTime()
                                  )
        inv_geonet = merge_inventories(inv_geonet, inv_hold)

    return inv_geonet[0]


def geonet_south_island(level="station"):
    """
    Return geonet broadband stations in an obspy network

    :type level: str
    :param level: level to propogate network creation
    :rtype: obspy.core.inventory.network.Network
    """
    c = Client("GEONET")

    # Extended North Island coverage, Kaikoura to East Cape
    lat_min = -48.
    lat_max = -40.
    lon_min = 165.
    lon_max = 175.

    starttime = UTCDateTime("2000-01-01")

    # GeoNet Broadband instruments with standard station code Z
    inv_geonet = c.get_stations(network="NZ", station="*Z", channel="HH?",
                                minlatitude=lat_min, maxlatitude=lat_max,
                                minlongitude=lon_min, maxlongitude=lon_max,
                                level=level, starttime=starttime,
                                endtime=UTCDateTime()
                                )

    return inv_geonet[0]


def hobitss(level="station"):
    """
    Return hobitss broadband OBS sensors (trillium compacts)

    :type level: str
    :param level: level to propogate network creation
    :rtype: obspy.core.inventory.network.Network
    """
    c = Client("IRIS")
    inv_hobitss = c.get_stations(network='YH', station="LOBS*", location='',
                                 channel="HH?", level=level
                                 )
    return inv_hobitss[0]


def sahke(network_code="X2", level="station", comp_list=["N", "E", "Z"]):
    """
    SAHKE transect broadband stations aren't fetchable through fdsn webservies,
    build a network object from the available information that was collected
    through the SAHKE final report provided by Martha Savage.

    I'm not sure if the CMG's are 30s or 60s instruments or how much that 
    actually matters

    Notes from GNS SAHKE Report:

        Instruments and dataloggers:
        CMG3ESP: T004, LTN6
        CMG40T: LE4, T007, T010, T014, T016, T018, T020
        Dataloggers: Reftek-130s 

        3-2-2010 (2010-061): LE4 sampling rate changed from 40Hz to 100Hz
        This isn't relevant for the data that I have...

    NRL variations:
        3ESP:
            Natural Period: "100 s - 50 Hz", "120 s - 50 Hz", "30 s - 50 Hz", 
                            "60 s - 50 Hz"
            Sensitivity: "1500", "2000", "20000"
        40T:
            Natural Period: "100s - 50Hz", "10s - 100Hz", "1s - 100Hz",
                            "20s - 50Hz", "2s - 100Hz", "30s - 100Hz",
                            "30s - 50 Hz", "40s - 100Hz", "5s - 100Hz",
                            "60s - 100Hz", "60s - 50Hz"
            Sensitivity:   "1600", "2000", "20000", "800"
        RT130S:
            Gain: "1", "32"



    :type network_code: str
    :param network_code: chosen two value code used for the network
    :type level: str
    :param level: level to propogate network creation
    :type comp_list: list of str
    :param comp_list: components to create channels for
    :rtype: obspy.core.inventory.network.Network
    :return: obspy Network object with information propogated to chosen level
    """
    # station, location, start, stop, lat, lon, instr type
    station_info = np.array([
        ["LE4",  "", "2010-136", "2010-331", -41.3579, 175.6919, "40t"],
        ["LTN6", "LT", "2010-193", "2010-349", -41.1033, 175.3238, "3esp"],
        ["T004", "", "2010-088", "2010-255", -41.3403, 175.6688, "3esp"],
        ["T007", "", "2010-041", "2010-123", -41.3041, 175.6513, "40t"],
        ["T010", "T0", "2010-135", "2010-348", -41.2520, 175.5825, "40t"],
        ["T014", "", "2010-034", "2010-350", -41.2075, 175.5063, "40t"],
        ["T016", "", "2010-088", "2010-322", -41.1893, 175.4737, "40t"],
        ["T018", "", "2010-055", "2010-349", -41.1715, 175.3850, "40t"],
        ["T020", "", "2010-089", "2010-261", -41.1251, 175.3497, "40t"]
    ])

    # For setting the network timing
    starttimes = station_info[:, 2]
    endtimes = station_info[:, 3]

    unique_starts = [UTCDateTime(str(_)) for _ in np.unique(starttimes)]
    unique_ends = [UTCDateTime(str(_)) for _ in np.unique(endtimes)]

    min_starttime = min(unique_starts)
    max_endtime = max(unique_ends)

    # Elevations are not known
    default_elevation = 0.0
    default_depth = 0.0
    default_site = Site(name="SAHKE")

    # Create response information
    if level == "channel":
        nrl = NRL()
        responses = {
            "40t": nrl.get_response(
                sensor_keys=["Guralp", "CMG-40T", "60s - 50Hz", "2000"],
                datalogger_keys=["REF TEK", "RT 130S & 130-SMHR", "1", "100"]),
            "3esp": nrl.get_response(
                sensor_keys=["Guralp", "CMG-3ESP", "60 s - 50 Hz", "2000"],
                datalogger_keys=["REF TEK", "RT 130S & 130-SMHR", "1", "100"]),
        }

    # Add stations to Station objects
    stations = []
    for stalist in station_info:
        # Parse the list to avoid confusion with indices
        code = stalist[0]
        location = stalist[1]
        start_date = UTCDateTime(stalist[2])
        end_date = UTCDateTime(stalist[3])
        latitude = stalist[4]
        longitude = stalist[5]

        # Create channel level objects if required
        if level == "channel":
            channels = []
            for comp in comp_list:
                cha = Channel(code=f"HH{comp}", location_code=location,
                              start_date=start_date, end_date=end_date,
                              latitude=latitude, longitude=longitude,
                              elevation=default_elevation, depth=default_depth,
                              azimuth=0.0, dip=-90.0, sample_rate=100
                              )
                cha.response = responses[stalist[-1]]
                channels.append(cha)
        else:
            channels = None

        # Create the Station object
        station = Station(code=code, start_date=start_date, end_date=end_date,
                          latitude=latitude, longitude=longitude,
                          elevation=default_elevation, site=default_site,
                          creation_date=UTCDateTime(), channels=channels
                          )

        stations.append(station)

    # Build the network and place a description
    network = Network(code=network_code, start_date=min_starttime,
                      end_date=max_endtime, description="SAHKE BBTRANSECT",
                      stations=stations
                      )

    return network


def bannister(network_code="ZX", level="station", comp_list=["N", "E", "Z"]):
    """
    Stephen Bannister deployed a broadband network in the northern north island.
    Data was provided via personal communication, this list of station locations
    was generated based on Stephen's email exchange with Yoshi.

    Most stations werent provided with an end date

    :type network_code: str
    :param network_code: chosen two value code used for the network
    :type level: str
    :param level: level to propogate network creation
    :type comp_list: list of str
    :param comp_list: components to create channels for
    :rtype: obspy.core.inventory.network.Network
    :return: obspy Network object with information propogated to chosen level
    """
    # station, lat, lon, depth, start, end
    station_info_zx = np.array([
        ['GA01', -39.0331, 177.8549, 33.000, 2011305, 9999001],
        ['GA02', -38.7793, 177.8695, 220.000, 2011305, 9999001],
        ['GA03', -38.6059, 177.9959, 46.000, 2011305, 9999001],
        ['GA04', -38.5529, 178.1548, 210.000, 2011305, 9999001],
        ['GA05', -38.5844, 177.8048, 85.000, 2011305, 9999001],
        ['GA06', -39.1110, 177.9161, 159.000, 2012001, 9999001],
        ['GA07', -38.4322, 178.0478, 569.000, 2012001, 9999001],
        ['GA08', -38.8312, 177.7095, 358.000, 2012001, 9999001],
        ['GA09', -38.5194, 177.9363, 107.000, 2013152, 9999001],
        ['GA10', -38.9087, 177.4627, 57.000, 2013152, 9999001],
        # ['GA11', -38.4932, 177.6711, 224.000, 2014060, 9999001],
        ])

    station_info_z8 = np.array([
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
        # We don't have any data for these stations
        # ['HD50', -38.3396, 176.2687, 361.000, 2015305, 9999001],
        # ['HD51', -38.2435, 176.2483, 380.000, 2015305, 9999001],
        # ['HD53', -38.2659, 176.4812, 362.000, 2015335, 9999001],
        # ['HD54', -38.2815, 176.3781, 480.000, 2015335, 9999001],
        # ['HD55', -38.2815, 176.5605, 469.000, 2015335, 9999001],
        # ['HD56', -38.3184, 176.5538, 439.000, 2015335, 9999001],
        # ['HD57', -38.1872, 176.3714, 520.000, 2015335, 9999001],
        # ['HD58', -38.3019, 176.3053, 492.000, 2015335, 9999001],
        # ['HD59', -38.2815, 176.3781, 440.000, 2015335, 9999001],
        # ['HD60', -38.2863, 176.1992, 389.000, 2016061, 9999001]
    ])

    # Pick which network to be exported
    if network_code == "ZX":
        station_info = station_info_zx
    elif network_code == "Z8":
        station_info = station_info_z8

    # Elevations are not known
    default_elevation = 0.0
    default_site = Site(name="BAN")

    # For setting the network timing
    unique_starts = [UTCDateTime(str(_)) for _ in np.unique(station_info[:, 4])]

    # Ignore the random value given for unknown endtimes
    unique_ends = [UTCDateTime(str(_)) for _ in
                   np.unique(station_info[:, 5])[:-1]]

    min_starttime = min(unique_starts)
    # If no actual endtimes, manual set the endtime
    try:
        max_endtime = max(unique_ends)
    except ValueError:
        max_endtime = UTCDateTime("2015-01-01")

    # Add stations to objects. Some endtimes are not specified
    stations = []
    for stalist in station_info:
        code = stalist[0]
        latitude = stalist[1]
        longitude = stalist[2]
        depth = stalist[3]
        start_date = UTCDateTime(str(stalist[4]))
        endtime = stalist[5]

        # Deal with endtimes, set to the latest available endtime in the list
        if endtime == 9999001:
            end_date = max(unique_ends)
        else:
            end_date = UTCDateTime(endtime)

        # Create channel level objects if required
        if level == "channel":
            channels = []
            for comp in comp_list:
                cha = Channel(code=f"HH{comp}", location_code="10",
                              start_date=start_date, end_date=end_date,
                              latitude=latitude, longitude=longitude,
                              elevation=default_elevation, depth=depth,
                              azimuth=0.0, dip=-90.0, sample_rate=100
                              )
                channels.append(cha)
        else:
            channels = None

        # Create the station object
        station = Station(code=code, latitude=latitude, longitude=longitude,
                          elevation=default_elevation, start_date=start_date,
                          end_date=end_date, site=default_site,
                          creation_date=UTCDateTime(), channels=channels
                          )
        stations.append(station)

    # Create the network object
    network = Network(code=network_code, start_date=min_starttime,
                      end_date=max_endtime,
                      description="Stephen Bannister's Broadband Network",
                      stations=stations
                      )

    return network


def beacon(network_code="2P", level="station", comp_list=["N", "E", "Z"]):
    """
    Create Beacon network data from scratch.
    Station information taken from the Site and Sensor field deployment notes
    kept on a shared Google Drive with Yoshi, Jonathan and myself.

    Updated Jun 23, 2022

    .. note::
        Start and end times are based on the MiniSEED data files

    :type network_code: str
    :param network_code: chosen two value code used for the network
    :type level: str
    :param level: level to propogate network creation
    :type comp_list: list of str
    :param comp_list: components to create channels for
    :rtype: obspy.core.inventory.network.Network
    :return: obspy Network object with information propogated to chosen level
    """
    # Station name, Abbreviation, Code, Lat, Lon, Start, End, Instrument type
    station_info = np.array([
        ["Pori Rd", "PORI", "RD01", "-40.55475083", "175.9710354",
            "2017-07-18T02:46:10.400000Z", 
            "2019-01-19T04:34:31.600000Z", "60s"],
        ["Angora Rd", "ANGR", "RD02", "-40.45974293", "176.4750588",
            "2017-07-18T00:06:19.090000Z", 
            "2019-01-19T13:36:37.970000Z", "60s"],
        ["Te Uri Rd", "TURI", "RD03", "-40.2656269", "176.3828498",
            "2017-07-18T00:10:35.140000Z", 
            "2019-02-22T05:19:15.270000Z", "30s"],
        ["Porangahau", "PORA", "RD04", "-40.2667317", "176.6344719",
            "2017-07-18T00:18:51.410000Z", 
            "2019-02-05T03:12:25.860000Z", "60s"],
        ["Manuhara Rd", "MNHR", "RD05", "-40.4689786", "176.2231874",
            "2017-07-18T04:08:06.500000Z", 
            "2019-02-22T02:45:06.830000Z", "30s"],
        ["Dannevirke", "DNVK", "RD06", "-40.2971794", "176.1663731",
            "2017-07-18T02:45:58.130000Z", 
            "2019-03-08T13:03:23.340000Z", "30s"],
        ["Waipawa", "WPAW", "RD07", "-39.9017124", "176.5370861",
            "2017-07-18T00:01:13.990000Z", 
            "2019-02-28T08:49:42.780000Z", "60s"],
        ["Raukawa", "RAKW", "RD08", "-39.7460611", "176.6205577",
            "2017-07-21T04:54:37.466100Z", 
            "2019-02-06T17:43:41.150000Z", "60s"],
        ["McNeill Hill", "MCNL", "RD09", "-39.4447675", "176.6974385",
            "2017-07-21T03:51:49.360000Z", 
            "2019-02-11T13:46:27.440000Z", "60s"],
        ["Cape Kidnappers", "CPKN", "RD10", "-39.64661592", "177.0765055",
            "2017-07-23T01:15:24.490000Z", 
            "2018-03-04T14:37:40.050000Z", "60s"],
        ["Kahuranaki", "KAHU", "RD11", "-39.78731589", "176.8624521",
            "2017-07-21T04:12:48.360000Z", 
            "2018-03-06T05:22:32.170000Z", "60s"],
        ["Kaweka Forest", "KWKA", "RD12", "-39.425214", "176.4228",
            "2017-07-21T04:22:08.830000Z", 
            "2019-02-04T18:04:08.470000Z", "30s"],
        ["Kereru", "KERE", "RD13", "-39.643259", "176.3768865",
            "2017-07-21T05:43:56.610000Z", 
            "2019-03-09T00:14:02.930000Z", "60s"],
        ["Pukenui", "PNUI", "RD14", "-39.9129963", "176.2001869",
            "2017-07-23T00:16:18.150000Z", 
            "2018-06-13T18:27:00.980000Z", "60s"],
        ["Waipukarau", "WPUK", "RD15", "-40.0627107", "176.4391311",
            "2017-07-23T00:10:44.120000Z", 
            "2019-02-11T00:09:50.690000Z", "60s"],
        ["Omakere", "OROA", "RD16", "-40.105341", "176.6804449",
            "2017-07-23T00:08:12.220000Z", 
            "2019-02-05T22:59:39.980000Z", "60s"],
        ["Te Apiti Rd", "TEAC", "RD17", "-39.90868978", "176.9561896",
            "2017-09-25T02:10:21.585100Z", 
            "2018-03-02T05:51:56.420000Z", "30s"],  # no sensor number, no instr type
        ["River Rd", "RANC", "RD18", "-39.929775", "176.7039773",
            "2017-09-25T04:55:05.610000Z", 
            "2019-01-18T20:12:29.850000Z", "30s"],
        ["Matapiro Rd", "MATT", "RD19", "-39.5796128", "176.6449024",
            "2018-03-13T20:31:38.610000Z", 
            "2018-06-23T16:11:59.100000Z", "30s"],  # same instr. as RD10
        ["Kahuranaki", "KAHU2", "RD20", "-39.79385769", "176.8758813",
            "2018-03-13T04:39:43.610000Z", 
            "2018-08-26T19:00:36.390000Z", "30s"],  # same instr. as RD11
        ["Te Apiti Rd", "TEAC2", "RD21", "-39.913152", "176.946881",
            "2018-03-10T09:22:33.610000Z", 
            "2019-01-22T22:57:49.770000Z", "30s"],  # same instr. as RD17
        ["Castlepoint", "CAPT", "RD22", "-40.910278", "176.199167",
            "2019-01-01T00:00:00.000000Z", 
            "2019-01-09T00:41:31.120000Z", "60s"],  # unknown sensor number
    ])

    # For setting the network timing
    starttimes = station_info[:, 5]
    endtimes = station_info[:, 6]

    unique_starts = [UTCDateTime(str(_)) for _ in np.unique(starttimes)]
    unique_ends = [UTCDateTime(str(_)) for _ in np.unique(endtimes)]

    min_starttime = min(unique_starts)
    max_endtime = max(unique_ends)

    # Elevations are not known
    default_elevation = 0.0
    default_depth = 0.0

    # Response is the same for all stations. Response information was provided
    # through personal correspondance to GeoNet site selection scientist
    # Jonathan Hanson, but could also be ascertained from the instrument type
    # and datalogger type
    if level == "channel":
        nrl = NRL()
        responses = {
            "30s": nrl.get_response(
                sensor_keys=["Guralp", "CMG-40T", "30s - 50 Hz", "800"],
                datalogger_keys=["Nanometrics", "Taurus", "16 Vpp (1)",
                                 "Low (default)", "Off", "100"]
            ),
            "60s": nrl.get_response(
                sensor_keys=["Guralp", "CMG-40T", "60s - 50Hz", "800"],
                datalogger_keys=["Nanometrics", "Taurus", "16 Vpp (1)",
                                 "Low (default)", "Off", "100"]
            )
        }
        sensors = {
            "30s": Equipment(type="sensor", manufacturer="Guralp",
                             model="CMG-40T", 
                             description=f"30s (Broadband) 30s-50Hz 800 V/m/s"),
            "60s": Equipment(type="sensor", manufacturer="Guralp",
                             model="CMG-40T",
                             description=f"60s (Broadband) 60s-50Hz 800 V/m/s"),
            }
        data_logger = Equipment(type="data_logger", manufacturer="Nanometrics",
                                model=f"Taurus",
                                description="16 Vpp (gain 1), 100 sps, low "
                                            "impedance, DC removal filter off")

    # Add stations to objects
    stations = []
    for stalist in station_info:
        # Parse the station information
        name = stalist[0]  # e.g. Castlepoint
        nickname = stalist[1]  # e.g. CAPT
        code = stalist[2]  # e.g. RD22
        latitude = float(stalist[3])
        longitude = float(stalist[4])
        start_date = UTCDateTime(stalist[5])
        end_date = UTCDateTime(stalist[6])
        sensor_type = stalist[7]

        # Create channel level objects if required
        if level == "channel":
            channels = []
            for comp in comp_list:
                cha = Channel(code=f"HH{comp}", location_code="10",
                              start_date=start_date, end_date=end_date,
                              latitude=latitude, longitude=longitude,
                              elevation=default_elevation, depth=default_depth,
                              azimuth=0.0, dip=-90.0, sample_rate=100,
                              sensor=sensors[sensor_type], 
                              data_logger=data_logger
                              )
                # Attach the response
                cha.response = responses[sensor_type]
                channels.append(cha)
        else:
            channels = None

        # Create the site object to provide information on the site location
        site = Site(name=nickname, description=name)

        # Create the station object
        station = Station(code=code, latitude=latitude, longitude=longitude,
                          elevation=default_elevation, start_date=start_date,
                          end_date=end_date, site=site,
                          creation_date=UTCDateTime(), channels=channels
                          )
        stations.append(station)

    # Create the network object
    network = Network(code=network_code, start_date=min_starttime,
                      end_date=max_endtime,
                      description="Broadband East Coast Network",
                      stations=stations
                      )
    return network


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
    write_to = "2P_BEACON_Dataless.xml"
    export_to_specfem = False
    export_to_seed_fmt = False
    plot = False

    # Create the Inventory
    master_inventory = Inventory(
        networks=[
            # geonet_south_island(level=level),
            # geonet_north_island(level=level),
            # hobitss(level=level),
            # sahke(level=level),
            # bannister(network_code="ZX", level=level),
            # bannister(network_code="Z8", level=level),
            beacon(level=level)
        ], source="PYATOA")
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
                              outfile="./master_inventory.png"
                              )
