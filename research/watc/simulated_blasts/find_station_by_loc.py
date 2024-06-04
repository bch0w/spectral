"""
Find a station based on a given lat/lon with a given wiggle room

.. rubric::

    $ python find_station_by_loc.py <lat> <lon>
"""
import sys
from obspy.clients.fdsn import Client
import pdb

lat_target = float(sys.argv[1])
lon_target = float(sys.argv[2])
try:
    wiggle_room = float(sys.argv[3])
except IndexError:
    wiggle_room = 0.25  # degrees

c = Client("IRIS")
inv = c.get_stations(channel="*", minlatitude=lat_target-wiggle_room,
                     maxlatitude=lat_target+wiggle_room, 
                     minlongitude=lon_target-wiggle_room, 
                     maxlongitude=lon_target+wiggle_room
                     )
print(inv)
pdb.set_trace()
