"""create stations for simulation runs with easy to read formatting
"""
import os
import glob

def true_if_outside_bounds(lat,lon,corners):
    """bool tell me if station falls outside the map bounds
    :type lat/lon: float
    :param lat/lon: object location
    :type corners: list of floats
    :param corners: map bounds [lat_bot,lat_top,lon_left,lon_right] 
                        e.g.   [-42,-36,173,178
    """
    lat_bot,lat_top,lon_left,lon_right = corners

    if (lon < lon_left) or (lon > lon_right):
        return True
    elif (lat < lat_bot) or (lat > lat_top):
        return True
    else:
        return False

def make_kupe_list():
    path = '/Users/chowbr/Documents/subduction/spectral/common/DATA/KUPEDATA/STATION_LISTS/ALL_STATIONS'
    with open(path,'r') as f:
        lines = f.readlines()
        
    linetemp = "{s:>8}{n:>6}{la:15.6f}{lo:15.6f}{el:8.2f}{de:8.2f}\n"
    corners = [-42.5,-36.957,173,179.499]
    with open(path+'_NEW','w') as f:
        for l in lines:
            l = l.strip().split()
            try:
                sta,net,lat,lon,_,dep,elv = l
            except ValueError:
                sta,net,lat,lon,dep,elv = l
            if net != 'NZ':
                continue
            lat = float(lat)
            lon = float(lon)
            dep = float(dep)
            elv = float(elv)
            if true_if_outside_bounds(lat,lon,corners):
                print(sta)
                continue
            f.write(linetemp.format(s=sta,n=net,la=lat,lo=lon,el=elv,de=dep))
            
def make_mapmaker_stations():
    from obspy.clients.fdsn import Client

    lat_lon=[-42.5007,-36.9488,172.9998,179.5077]
    c = Client("GEONET")
    inv = c.get_stations(network='NZ',
                        station='*Z',
                        channel='?H?',
                        minlatitude=lat_lon[0],
                        maxlatitude=lat_lon[1],
                        minlongitude=lat_lon[2],
                        maxlongitude=lat_lon[3],
                        level="station")
    