"""
Quick export script to get seismic stations as KML files for Google Earth
"""
from obspy import Inventory
from obspy.clients.fdsn import Client

c = Client("GEONET")
networks = ["NZ"]
inv = Inventory()
for network in networks:
    inv += c.get_stations(network=network, station="*", location="*", 
                          channel="?H?", level="station")

print(inv)
inv.write(f"stations_{'_'.join(networks)}.kml", format="KML")
