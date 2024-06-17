"""
Given station names, find the exact locations and create a StationXML file
"""
from obspy import Inventory, UTCDateTime
from obspy.clients.fdsn import Client
from pyatoa.utils.srcrcv import merge_inventories


c = Client("IRIS")
inv = Inventory()
err = []
start = UTCDateTime("2024-01-01T00:00:00")  # to ensure we get active stations

with open("data/ims_primary_seismic.csv") as f:
# with open("data/ims_auxiliary_seismic.csv") as f:
    lines = f.readlines()
    for line in lines[1:]:
        print(line.split(",")[0], end=" ")
        codes = line.split(",")[4:]
        codes = [_.strip().replace('"', "") for _ in codes]
        check = 0
        for code in codes:
            if not code or code == "":
                continue
            net, sta = code.split(".")
            inv_ = None
            try:
                inv_ = c.get_stations(network=net, station=sta, 
                                      channel="HH?,BH?,EH?,SH?",
                                      level="channel",
                                      starttime=start)
                inv = merge_inventories(inv, inv_)
                check += 1
            except Exception as e:
                continue
        if not check:
            print("ERROR")
        print("")

import pdb;pdb.set_trace()




