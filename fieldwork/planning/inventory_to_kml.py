"""
Generate Google Earth KML files from ObsPy Inventory objects
"""
from obspy.clients.fdsn import Client

c = Client("IRIS")
# Bounding box on Alaska, getting seismometers
inv = c.get_stations(channel="?H?", minlatitude=50.8, maxlatitude=72.1, 
                     minlongitude=-179., maxlongitude=-150, level="station")
inv.write("inventory.kml", format="KML")

