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


def geonet_north_island(level="station"):
    """
    Return geonet broadband stations in an obspy network

    :rtype: obspy.core.inventory.network.Network
    """
    c = Client("GEONET")

    # Extended North Island coverage, Kaikoura to East Cape
    lat_min = -42.8
    lat_max = -37.1
    lon_min = 173
    lon_max = 179.4

    starttime = UTCDateTime('2000-01-01')

    # GeoNet Broadband instruments
    inv_geonet = c.get_stations(network='NZ', station='*Z', channel='HH?',
                                minlatitude=lat_min, maxlatitude=lat_max,
                                minlongitude=lon_min, maxlongitude=lon_max,
                                level=level, starttime=starttime,
                                endtime=UTCDateTime()
                                )

    return inv_geonet[0]


def hobitss(level="station"):
    """
    Return hobitss broadband OBS sensors (trillium compacts)

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
    via personal communication at VUW

    NO RESPONSE CREATED

    :rtype: obspy.core.inventory.network.Network
    """
    # station, start, stop, lat, lon
    station_info = np.array([
        ["LE4", "2010-136", "2010-331", -41.3579, 175.6919],
        ["LTN6", "2010-193", "2010-349", -41.1033, 175.3238],
        ["T004", "2010-088", "2010-255", -41.3403, 175.6688],
        ["T007", "2010-041", "2010-123", -41.3041, 175.6513],
        ["T010", "2010-135", "2010-348", -41.2520, 175.5825],
        ["T014", "2010-034", "2010-350", -41.2075, 175.5063],
        ["T016", "2010-088", "2010-322", -41.1893, 175.4737],
        ["T018", "2010-055", "2010-349", -41.1715, 175.3850],
        ["T020", "2010-089", "2010-261", -41.1251, 175.3497]
    ])

    # For setting the network timing
    starttimes = station_info[:, 1]
    endtimes = station_info[:, 2]

    unique_starts = [UTCDateTime(str(_)) for _ in np.unique(starttimes)]
    unique_ends = [UTCDateTime(str(_)) for _ in np.unique(endtimes)]

    min_starttime = min(unique_starts)
    max_endtime = max(unique_ends)

    # Elevations are not known
    default_elevation = 0.0
    default_depth = 0.0
    default_site = Site(name="SAHKE")

    # Add stations to Station objects
    stations = []
    for stalist in station_info:
        # Parse the list to avoid confusion with indices
        code = stalist[0]
        start_date = UTCDateTime(stalist[1])
        end_date = UTCDateTime(stalist[2])
        latitude = stalist[3]
        longitude = stalist[4]

        # Create channel level objects if required
        if level == "channel":
            channels = []
            for comp in comp_list:
                cha = Channel(code=f"HH{comp}", location_code="",
                              start_date=start_date, end_date=end_date,
                              latitude=latitude, longitude=longitude,
                              elevation=default_elevation, depth=default_depth,
                              azimuth=0.0, dip=-90.0, sample_rate=100
                              )
                channels.append(cha)
        else:
            channels = None

        # Create the Station object
        station = Station(code=code, start_date=start_date, end_date=end_date,
                          latitude=latitude, longitude=longitude,
                          elevation=default_elevation, site=default_site,
                          creation_date=UTCDateTime(),channels=channels
                          )

        stations.append(station)

    # Build the network and place a description
    network = Network(code=network_code, start_date=min_starttime,
                      end_date=max_endtime,
                      description="SAHKE2 Broadband Transect",
                      stations=stations
                      )

    return network


def bannister(network_code="X1", level="station", comp_list=["N", "E", "Z"]):
    """
    Stephen Bannister deployed a broadband network in the northern north island.
    Data was provided via personal communication, this list of station locations
    was generated based on Stephen's email exchange with Yoshi.

    :return:
    """
    # station, lat, lon, depth, start, end
    station_info = np.array([
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
    ])

    # Elevations are not known
    default_elevation = 0.0
    default_site = Site(name="BAN")

    # For setting the network timing
    unique_starts = [UTCDateTime(str(_)) for _ in np.unique(station_info[:, 4])]

    # Ignore the random value given for unknown endtimes
    unique_ends = [UTCDateTime(str(_)) for _ in
                   np.unique(station_info[:, 5])[:-1]]

    min_starttime = min(unique_starts)
    max_endtime = max(unique_ends)

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
                cha = Channel(code=f"HH{comp}", location_code="",
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


def beacon(network_code="XX", level="station", comp_list=["N", "E", "Z"]):
    """
    Create Beacon network data from scratch.
    Station information taken from the Site and Sensor field deployment notes
    kept on a shared Google Drive with Yoshi, Jonathan and myself.

    :return:
    """
    # Station name, Abbreviation, Code, Lat, Lon, Start, End
    station_info = np.array([
        ['Pori Rd', 'PORI', 'RD01', '-40.55475083', '175.9710354',
         '2017-07-19', '2019-04-04'],
        ['Angora Rd', 'ANGR', 'RD02', '-40.45974293', '176.4750588',
         '2017-07-19', '2019-04-04'],
        ['Te Uri Rd', 'TURI', 'RD03', '-40.2656269', '176.3828498',
         '2017-07-20', '2019-04-04'],
        ['Porangahau', 'PORA', 'RD04', '-40.2667317', '176.6344719',
         '2017-07-20', '2019-04-04'],
        ['Manuhara Rd', 'MNHR', 'RD05', '-40.4689786', '176.2231874',
         '2017-07-20', '2019-04-05'],
        ['Dannevirke', 'DNVK', 'RD06', '-40.2971794', '176.1663731',
         '2017-07-24', '2019-04-02'],
        ['Waipawa', 'WPAW', 'RD07', '-39.9017124', '176.5370861',
         '2017-07-24', '2019-04-02'],
        ['Raukawa', 'RAKW', 'RD08', '-39.7460611', '176.6205577',
         '2017-07-24', '2019-04-02'],
        ['McNeill Hill', 'MCNL', 'RD09', '-39.4447675', '176.6974385',
         '2017-07-25', '2019-04-03'],
        ['Cape Kidnappers', 'CPKN', 'RD10', '-39.64661592', '177.0765055',
         '2017-07-25', '2018-03-13'],
        ['Kahuranaki', 'KAHU', 'RD11', '-39.78731589', '176.8624521',
         '2017-07-25', '2018-03-13'],
        ['Kaweka Forest', 'KWKA', 'RD12', '-39.425214', '176.4228',
         '2017-07-26', '2019-05-03'],
        ['Kereru', 'KERE', 'RD13', '-39.643259', '176.3768865',
         '2017-07-26', '2019-04-03'],
        ['Pukenui', 'PNUI', 'RD14', '-39.9129963', '176.2001869',
         '2017-07-26', '2018-09-08'],
        ['Waipukarau', 'WPUK', 'RD15', '-40.0627107', '176.4391311',
         '2017-07-27', '2019-04-02'],
        ['Omakere', 'OROA', 'RD16', '-40.105341', '176.6804449',
         '2017-07-27', '2019-04-04'],
        ['Te Apiti Rd', 'TEAC', 'RD17', '-39.90868978', '176.9561896',
         '2017-09-25', '2018-03-14'],
        ['River Rd', 'RANC', 'RD18', '-39.929775', '176.7039773',
         '2017-09-25', '2019-04-03'],
        ['Matapiro Rd', 'MATT', 'RD19', '-39.5796128', '176.6449024',
         '2018-03-14', '2018-06-25'],
        ['Kahuranaki', 'KAHU2', 'RD20', '-39.79385769', '176.8758813',
         '2018-03-13', '2018-09-03'],
        ['Te Apiti Rd', 'TEAC2', 'RD21', '-39.913152', '176.946881',
         '2018-03-14', '2019-04-03'],
        ['Castlepoint', 'CAPT', 'RD22', '-40.910278', '176.199167',
         '2018-07-20', '2019-05-05']
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
        response = nrl.get_response(
            sensor_keys=['Guralp', 'CMG-40T', '60s - 50Hz', '800'],
            datalogger_keys=['Nanometrics', 'Taurus', '40 Vpp (0.4)',
                             'Low (default)', '1 mHz', '100']
        )

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

        # Create channel level objects if required
        if level == "channel":
            channels = []
            for comp in comp_list:
                cha = Channel(code=f"HH{comp}", location_code="",
                              start_date=start_date, end_date=end_date,
                              latitude=latitude, longitude=longitude,
                              elevation=default_elevation, depth=default_depth,
                              azimuth=0.0, dip=-90.0, sample_rate=100
                              )
                # Attach the response
                cha.response = response
                channels.append(cha)
        else:
            channels = None

        # Create the site object to provide information on the site location
        site = Site(name=nickname, description=name)

        # Create the station object
        station = Station(code=code, latitude=latitude,
                          longitude=longitude, elevation=default_elevation,
                          start_date=start_date, end_date=end_date, site=site,
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


def generate_master_list():
    """
    Call all functions in this script to create a master list of station data

    :return:
    """
    level = "channel"
    geonet_inv = geonet_north_island(level=level)
    hobitss_inv = hobitss(level=level)
    sahke_inv = sahke(level=level)
    bannister_inv = bannister(level=level)
    beacon_inv = beacon(level=level)

    return Inventory(networks=[geonet_inv, hobitss_inv, sahke_inv,
                               bannister_inv, beacon_inv],
                     source="PYATOA"
                     )


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


def export_pyatoa(inv):
    """
    Pyatoa requires a "seed" convention directory structure containing
    individual RESPONSE xml files. Create that here based on Obspy inventory

    :param inventory:
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
    master_inventory = generate_master_list()
    master_inventory.plot(projection='local', resolution='l',
                          continent_fill_color='w', water_fill_color='w',
                          label=True, color_per_network=True,
                          outfile="./master_inventory.png"
                          )
    master_inventory.write('./master_inventory.xml', format='STATIONXML')
    export_pyatoa(master_inventory)
