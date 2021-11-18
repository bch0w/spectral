"""
Write a catalog into a .csv file that can be read by Carl's seis_decluster
"""
import sys
from obspy import read_events

cat = read_events(sys.argv[1])
with open(f"{sys.argv[1]}.csv", "w") as f:
    f.write("event_id,time,magnitude,depth_km,latitude,longitude\n")
    for event in cat:
        eid = event.resource_id.id.split("/")[1]
        print(eid)
        tme = event.preferred_origin().time
        lat = event.preferred_origin().latitude
        lon = event.preferred_origin().longitude
        dpt = event.preferred_origin().depth * 1E-3
        mag = event.preferred_magnitude().mag
        f.write(f"{eid},{tme},{mag},{dpt},{lat},{lon}\n")


