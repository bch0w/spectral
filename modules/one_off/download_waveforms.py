from __future__ import print_function
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# set parameters here
network = "NZ"
location = "*"
channel = "B??"  # or "HHZ", "??E", "*N" etc.
station = "????"
output = "VEL"  # also "DISP" or "ACC"
starttime = UTCDateTime("2016-11-22T00:19:43.00Z")
endtime = starttime + 300
rotate_to_radial_transverse = False
c = Client("GEONET")

# retrieve all station codes dynamically, e.g. 4 letter strong motion stations
station_inv = c.get_stations(network=network, station=station,
                             location=location, channel=channel,
                             starttime=starttime,
                             endtime=endtime, level="station")
station_list = []
for station_ in station_inv[0]:
    station_list.append(station_.code)

# OR manually set station codes in a list here
# station_list = ["PGFS"]  # or ["BKZ", "PUZ", "KNZ"] etc.


# download data for each station separately
for station in station_list:
    print(station, end="...")
    st = c.get_waveforms(network=network, station=station, location=location,
                         channel=channel, starttime=starttime, endtime=endtime,
                         attach_response=True)
    st.remove_response(output=output)

    # trim because get_waveforms() sometimes returns non matching lengths
    st.trim(starttime=starttime, endtime=endtime)

    # check if components are in 1,2,Z, if so rotate by azimuth and inclination
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
        tr.write(filename=filename, format="SAC")
    print("saved")

