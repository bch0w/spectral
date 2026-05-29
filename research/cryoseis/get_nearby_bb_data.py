import sys
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# PARAMETERS
network = "AK"
station = "PAX"  # PS10
location ="*"
channel = "BH?"

save_to = "/Users/chow/Work/research/gulkanaseis24/akdata"
    
c = Client("IRIS")
for i in range(251, 259):
    starttime = UTCDateTime(f"2024-{i}T00:00:00")
    endtime = UTCDateTime(starttime) + (24 * 60 * 60 - .01)

    st = c.get_waveforms(network=network, station=station, location=location,
                         channel=channel, starttime=starttime, endtime=endtime, 
                         attach_response=True)
    print(st)
    st.remove_response()
    for tr in st:
        _time = tr.stats.starttime
        tr.write(f"{save_to}/{tr.get_id()}.{_time.year}.{_time.julday}", 
                 format="MSEED")

