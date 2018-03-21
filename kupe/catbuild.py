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

def events_by_radius():
    c = Client("GEONET")
    new_zealand = [-50,-35,165,180]
    lat_lon = new_zealand
    # all M5-6 events in New Zealand
    cat = c.get_events(starttime="2006-01-01T00:00:00",
                        endtime=UTCDateTime(),
                        minmagnitude=5,
                        maxmagnitude=6,
                        minlatitude=lat_lon[0],
                        maxlatitude=lat_lon[1],
                        minlongitude=lat_lon[2],
                        maxlongitude=lat_lon[3],
                        orderby="magnitude")
    for net in inv:
        for sta in net:
            lat = sta.latitude
            lon = sta.longitude
            start = sta.creation_date
            end = sta.termination_date
            if not end:
                end = UTCDateTime()
            # 4 <= M <= 5
            cat += c.get_events(starttime=start,
                                endtime=end,
                                minmagnitude=4,
                                maxmagnitude=5,
                                latitude=lat,
                                longitude=lon,
                                maxradius=3)
            # 3 <= M <= 4
            cat += c.get_events(starttime=start,
                                endtime=end,
                                minmagnitude=3,
                                maxmagnitude=4,
                                latitude=lat,
                                longitude=lon,
                                maxradius=2)
            # 2 <= M <= 3
            cat += c.get_events(starttime=start,
                                endtime=end,
                                minmagnitude=2,
                                maxmagnitude=3,
                                latitude=lat,
                                longitude=lon,
                                maxradius=1)
