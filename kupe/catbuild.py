"""to build up a catalog of earthquakes useable for seismic tomography

potential workflow outline
-collect all stations on the north island, as well as temporary deployments,
offshore deployments, rdf etc.
-for each station, collect stationxml information into a single xml file, and
numpy arrays containing lat/lon of each station
-collect all m4 or m5 events in the area; function ability to manually look at
waveforms and assess quality
-collect all <m3 events for a radius surrounding each station for more events,
using lat/lon as inputs to event getter
-collect all quakeml information into single quakeml file, somehow organize
which events belong to which stations (excel file?)
-run detection algorithm on all station data?
-run flexwin on all picks for optimum window choice?
"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

def one_time_collect_stations():
    # geonet permanent stations
    c = Client("GEONET")
    north_island = [-42,-34,173,180]
    lat_lon = north_island
    inv = c.get_stations(network='NZ',
                        station='*',
                        channel='HH*',
                        minlatitude=lat_lon[0],
                        maxlatitude=lat_lon[1],
                        minlongitude=lat_lon[2],
                        maxlongitude=lat_lon[3],
                        level="response")
    # hobitss stations
    c = Client("IRIS")
    inv += c.get_stations(network='YH',
                        station="LOBS*",
                        location='',
                        level="response")
    inv += c.get_stations(network='YH',
                        station="EBS*",
                        location='',
                        level="response")

    inv.write("station_massive.xml",format="STATIONXML")

def collect_events():
    c = Client("GEONET")
    new_zealand = [-50,-35,165,180]
    kaikoura_to_east_cape = [-43,-37,172,180]
    lat_lon = kaikoura_to_east_cape
    cat = c.get_events(starttime="2005-01-01T00:00:00",
                        endtime=UTCDateTime(),
                        minmagnitude=4.5,
                        maxmagnitude=6,
                        maxdepth=80,
                        minlatitude=lat_lon[0],
                        maxlatitude=lat_lon[1],
                        minlongitude=lat_lon[2],
                        maxlongitude=lat_lon[3],
                        orderby="magnitude")
