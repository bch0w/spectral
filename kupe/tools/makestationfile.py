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
    
    