"""
Generate Google Earth KML files from ObsPy Inventory objects
"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

kwargs = dict(
        channel="?H?", 
        starttime="2010-01-01T00:00:00",
        minlatitude=50.8, 
        maxlatitude=72.1,
        minlongitude=-179., 
        maxlongitude=-140, 
        level="channel"
        )

c = Client("IRIS")

# Bounding box on Alaska, getting seismometers
inv = c.get_stations(**kwargs)
inv.write("Inventory.kml", format="KML")

