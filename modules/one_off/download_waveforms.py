from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# set parameters here
network = "NZ"
station = "BKZ"
location = "*"
channel = "HH?"  # or HHZ, ??E, *N etc.
output = "VEL"  # also "DISP" or "ACC"
starttime = UTCDateTime("2018-01-01T00:00:00")
endtime = starttime + 60*60  # or UTCDateTime("2018-01-01T01:00:00")

# data downloaded here
c = Client("GEONET")
st = c.get_waveforms(network=network, station=station, location=location,
                     channel=channel, starttime=starttime, endtime=endtime,
                     attach_response=True)
st.remove_response(output=output)
for tr in st:
    # filename will be saved e.g. NZ.BFZ.10.HHZ.2018.001.SAC, can also be
    # set manually
    filename="{id}.{year}.{jday:0>3}.SAC".format(
        id=tr.get_id(), year=tr.stats.starttime.year,
        jday=tr.stats.starttime.julday)
    tr.write(filename=filename, format="SAC")
