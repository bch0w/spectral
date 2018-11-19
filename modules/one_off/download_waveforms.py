"""
script to download seismic waveforms from the FDSN webservice via ObsPy.
Functionality available to search internal GeoNet directories for data

Parameters are set in script in the denoted section. Currently the script is
written with the intent to download data from GeoNet, but ObsPy allows
data download from any databases accessible via FDSN.

Misc. Information:
-For convenience, output written directly into the current working directory
-Files are saved individually by component
-Filenames generated follow the template:
    network.station.location.channel.year.julianday.format
-Data given in 1/2/Z components (e.g. accelerometers) are automatically rotated
 to an E/N/Z coordinate system

Relevant Obspy Documentation:
+ For available data formats:
    https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.write.html
+ For information regarding FDSN clients:
    https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.get_waveforms.html
+ For a list of available clients:
    https://docs.obspy.org/packages/obspy.clients.fdsn.html#module-obspy.clients.fdsn
"""
from __future__ import print_function
import os
import glob
from obspy import read, UTCDateTime
from obspy.core.stream import Stream
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.client import FDSNNoDataException


def overlapping_days(starttime, endtime):
    """
    stolen and edited from Pyaota

    Helper function to return a list of julian days based on a given
    origin time with a specific padding on either side. used to catch if an
    origin time sits too close to midnight and two days need to be fetched
    """
    if starttime.julday != endtime.julday:
        return [starttime.julday, endtime.julday]
    else:
        return [starttime.julday]


def fetch_by_directory(network, station, location, channel, starttime, endtime):
    """
    stolen and edited from Pyaota

    waveform directory structure formatting is hardcoded in here, it is
    assumed that data is saved as miniseeds in the following structure:

    path/to/data/{YEAR}/{NETWORK}/{STATION}/{CHANNEL}*/{FID}
    """
    # HARDCODED GEONET LOCAL DIRECTORY PATH
    path_ = '/geonet/seismic/'

    dir_structure = '{year}/{net}/{sta}/{cha}*'
    file_template = '{net}.{sta}.{loc}.{cha}*{year}.{jday:0>3}'
    jdays = overlapping_days(starttime=starttime, endtime=endtime)

    if not os.path.exists(path_):
        print("path '/geonet/seismic' does not exist")
        return None
    full_path = os.path.join(path_, dir_structure, file_template)
    pathlist = []
    for jday in jdays:
        pathlist.append(full_path.format(net=network, sta=station, cha=channel,
                                         loc=location, jday=jday,
                                         year=starttime.year)
                        )
    st = Stream()
    for fid in pathlist:
        for filepath in glob.glob(fid):
            st += read(filepath)
    if len(st) > 0:  # is this necessary?
        st.merge()
        st.trim(starttime=starttime, endtime=endtime)
        return st
    else:
        raise FileNotFoundError(
            "No waveforms found for {} for given directories".format(
                file_template.format(net=network, sta=station, cha=channel,
                                     loc=location, jday=jday,
                                     year=starttime.year))
        )


"""
======================== vv SET PARAMETERS HERE vv ===========================
:type network: str
:param network: network code, e.g. "NZ"
:type location: str
:param location: location code, e.g. "00", can also be wildcard, e.g. "*"
:type channel:  str
:param channel: channel code, e.g. "HHZ", wildcard okay, e.g. "??E", "*N", etc.
:type station: str
:param station: station code, e.g. "PUZ", wildcard okay
:type output: str
:param output: unit output, available "DISP", "VEL", "ACC" or None for raw data
:type output_format: str
:param output_format: file format to save waveforms, e.g. "SAC", "MSEED", etc.
:type starttime: UTCDateTime(str)
:param starttime: start of waveforms, must  be wrapped in UTCDateTime()
:type endtime: UTCDateTime(str) or UTCDateTime
:param endtime: end of waveforms, can be the the same format as starttime, or
    integer addition to starttime, i.e. endtime = starttime + 300 
:type rotate_to_radial_transverse: bool
:param rotate_to_radial_transverse: leave in NEZ (or 12Z) or rotate to RTZ
:type c: obspy.clients.fdsn
:param c: client to download data from, e.g. GEONET, IRIS, GFZ
:type suppress_fdsn: bool
:param suppress_fdsn: skip searching via FDSN webservices, and only query 
    local geonet directories for data. Potentially faster if you know that data
    is not available on the FDSN webservices. (I'm not sure if FDSN or 
    internal fetching is faster, because internal fetching requires reading
    the entire day-length mseed file, sometimes two if starttimes near 00:00)
:type station_select_type: int
:param station_select_type: choose 1 for wildcard selection, 2 for manual list
:type station_list: list
:param station_list: list of strings containing station names, only relevatn if
    station_select_type is set to 2
"""

# station information
network = "NZ"
location = "*"
channel = "B??"
station = "????"

# event information
starttime = UTCDateTime("2018-10-30T02:13:41.649Z")
endtime = starttime + 300

# output information
output = "VEL"
rotate_to_radial_transverse = False
output_format = "SAC"

# FDSN Downloading
suppress_fdsn = True
c = Client("GEONET")

# ------------------------- vv CHOOSE STATIONS BY vv ---------------------------
station_select_type = 1 # or 2
# 1) Wildcard selection: will dynamically create a list of stations based on
#    the given parameters above.
if station_select_type == 1:
    station_inv = c.get_stations(network=network, station=station,
                                 location=location, channel=channel,
                                 starttime=starttime,
                                 endtime=endtime, level="station")
    station_list = []
    for station_ in station_inv[0]:
        station_list.append(station_.code)

# OR
# 2) Manual: set station codes in a list, ignores 'station' parameter above
if station_select_type == 2:
    station_list = ["ADCS"]  # or ["BKZ", "PUZ", "KNZ"] etc.

# ======================== ^^ SET PARAMETERS HERE ^^ ===========================

# Data download on a per-station basis begins here
print("Downloading data for {} station(s):".format(len(station_list)))
for i, station in enumerate(station_list):
    try:
        print("{}/{} {}".format(i + 1, len(station_list), station), end="...")
        # first try downloading via FDSN webservices, except no data available
        try:
            if suppress_fdsn:
                raise FDSNNoDataException(True)  # cheeky way to skip external
            st = c.get_waveforms(network=network, station=station,
                                 location=location, channel=channel,
                                 starttime=starttime, endtime=endtime,
                                 attach_response=True)
            st.trim(starttime=starttime, endtime=endtime)
            if output is not None:
                st.remove_response(output=output)
            print("via fdsn", end="...")

        # if failed or skipped external search, search internal GeoNet pathways
        except FDSNNoDataException:
            st = fetch_by_directory(network=network, station=station,
                                    location=location, channel=channel,
                                    starttime=starttime, endtime=endtime,
                                    )
            inv = c.get_stations(network=network, station=station,
                                 location=location, channel=channel,
                                 starttime=starttime,
                                 endtime=endtime, level="response")
            st.remove_response(inventory=inv, output=output)
            print("via local", end="...")

        # check if components are in 1,2,Z, if so rotate by azi. and inclin.
        # values retrieved in get_stations()
        st.sort()
        if len(st) == 3 and st[0].get_id()[-1] == "1":
            print("rotating 1/2/Z to N/E/Z", end="...")
            inv = c.get_stations(network=network, station=station,
                                 location=location, channel=channel,
                                 starttime=starttime, endtime=endtime,
                                 level="channel")
            st.rotate(method="->ZNE", inventory=inv)

        # if asked for, rotate from NEZ to RTZ
        if rotate_to_radial_transverse:
            print("rotating N/E/Z to R/T/Z", end="...")
            st.rotate(method="NE->RT")

        # write out streams
        for tr in st:
            # filename will be saved e.g. NZ.BFZ.10.HHZ.2018.001.SAC
            filename = "{id}.{year}.{jday:0>3}.SAC".format(
                id=tr.get_id(), year=tr.stats.starttime.year,
                jday=tr.stats.starttime.julday)
            tr.write(filename=filename, format=output_format)
        print("saved")
    except Exception as e:
        print("failed, for reason:\n\t{}".format(e))
        continue
