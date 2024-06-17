"""
Find a station based on a given coordinates from a CSV file
"""
import sys
import numpy as np
from obspy.clients.fdsn import Client

primary = np.loadtxt("ims_primary_seismic.csv", skiprows=1, usecols=[0,1,2,3],
                     delimiter=",", dtype="object")

auxiliary = np.loadtxt("ims_auxiliary_seismic.csv", skiprows=1,
                       usecols=[0,1,2,3], delimiter=",", dtype="object")


c = Client("IRIS")
output = {}
# for p in primary:
for p in auxiliary:
    inv = None
    wiggle_room = 0  # degrees
    sta_target, lat_target, lon_target, _ = p
    lat_target = float(lat_target)
    lon_target = float(lon_target)
    print(sta_target)
    output[sta_target] = []
    while True:
        try:
            inv = c.get_stations(channel="BH?,HH?", 
                                 minlatitude=lat_target-wiggle_room,
                                 maxlatitude=lat_target+wiggle_room, 
                                 minlongitude=lon_target-wiggle_room, 
                                 maxlongitude=lon_target+wiggle_room,
                                 level="station"
                                 )
            break
        except Exception as e:
            wiggle_room += .1
            if wiggle_room >= 2.5:
                break
    if inv:
        for net in inv:
            if net.code in ["SY"]:
                continue
            for sta in net:
                netsta = f"{net.code}.{sta.code}"
                if netsta not in output[sta_target]:
                    output[sta_target].append(f"{net.code}.{sta.code}")

for key, val in output.items():
    if not val:
        print(" ")
        continue
    elif val:
        print(", ".join(val))
