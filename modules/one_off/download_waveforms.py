"""
Short script to download seismic waveforms from the FDSN webservice via ObsPy.

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
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# ======================== vv SET PARAMETERS HERE vv ===========================
network = "NZ"
location = "*"
channel = "B??"  # or "HHZ", "??E", "*N" etc.
station = "????"
output = "VEL"  # "DISP", "VEL", "ACC" or None if you want raw data
starttime = UTCDateTime("2016-11-22T00:19:43.00Z")
endtime = starttime + 300
rotate_to_radial_transverse = False
output_format = "SAC"  # also e.g. "MSEED"
c = Client("GEONET")

# ------------------------- vv CHOOSE STATIONS BY vv ---------------------------
station_select_type = 1  # or 2
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
        print("{}/{} {}".format(i+1, len(station_list), station), end="...")
        st = c.get_waveforms(network=network, station=station,
                             location=location, channel=channel,
                             starttime=starttime, endtime=endtime,
                             attach_response=True)
        if output is not None:
            st.remove_response(output=output)

        # trim because get_waveforms() sometimes returns non matching lengths
        st.trim(starttime=starttime, endtime=endtime)

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
        print("failed, for reason:\n{}".format(e))
        continue
