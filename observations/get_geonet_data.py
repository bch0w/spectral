"""
A short Python script used to download GeoNet data
This requires ObsPy (obspy.org)
"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# vvv  SET PARAMETERS HERE vvv
network = "NZ"  # this is usually static
station = "BFZ"   
location = "*"  # GeoNet broadband locations are usually '10'
channel = "HH?"  # wildcard component to get 3 components
starttime = UTCDateTime("2012-01-01T00:00:00")  # Set to earthquake origintime
endtime = starttime + 300  # end of waveform, X seconds are starttime
remove_response = True  # if you want to remove the instrument response
# ^^^ SET PARAMETERS HERE ^^^

c = Client("GEONET")
st = c.get_waveforms(network=network, station=station, location=location,
                     channel=channel, starttime=starttime, endtime=endtime,
                     attach_response=True)
if remove_response:
    st.remove_response()



