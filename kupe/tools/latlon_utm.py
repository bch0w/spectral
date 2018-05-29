"""one off script for converting lat lon coordinates to UTM-60 projection
"""

def latlon_utm(lon_or_x,lat_or_y,zone=60,inverse=False):
    """convert latitude and longitude coordinates to UTM projection
    if inverse == False, latlon => UTM
    """
    from pyproj import Proj

    myProj = Proj(proj='utm',zone=zone,ellps='WGS84')
    x_or_lon,y_or_lat = myProj(lon_or_x,lat_or_y,inverse=inverse)
    
    return x_or_lon, y_or_lat
    
if __name__ == "__main__":
    LLC_lonlat = [173.3795,-42.4198]
    URC_lonlat = [178.8848,-37.4556]
    LLC_xy= latlon_utm(LLC_lonlat[0],LLC_lonlat[1])
    URC_xy= latlon_utm(URC_lonlat[0],URC_lonlat[1])
    
    
    
