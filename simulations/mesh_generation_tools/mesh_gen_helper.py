"""helper functions designed to facilitate the process of meshing and 
running simulations with kupe. separated into specific functionalities

determine_nxny: convet WGS84 to UTM and find optimum mesh size
lonlat_utm: convert WGS84 to UTM
"""

# ================================ MESHING ====================================
def lonlat_utm(lon_or_x,lat_or_y,zone=60,inverse=False):
    """convert latitude and longitude coordinates to UTM projection
    
    :type lon_or_x: float or int
    :param lon_or_x: longitude value in WGS84 or X in UTM-'zone' projection
    :type lat_or_y: float or int
    :param lat_or_y: latude value in WGS84 or Y in UTM-'zone' projection
    :type zone: int
    :param zone: UTM zone for conversion from WGS84
    :type inverse: bool
    :param inverse: if inverse == False, latlon => UTM, vice versa.
    """
    from pyproj import Proj
    projstr = ("+proj=utm +zone=60, +south +ellps=WGS84 +datum=WGS84"
                                                        " +units=m +no_defs")
    myProj = Proj(projstr)
    x_or_lon,y_or_lat = myProj(lon_or_x,lat_or_y,inverse=inverse)
    
    return x_or_lon, y_or_lat
    
    
def determine_nxny(LLC_lonlat,URC_lonlat):
    """given coordinates in WGS84, convert to UTM60 and return optimal nx and ny
    divisions for mesh generation. Each input variable should be a tuple
    gives a call to lonlat_utm
    :type LLC_lonlat: tuple
    :param LLC_lonlat: lower left hand corner [longitude, latitude]
    :type URC_lonlat: tuple
    :param URC_lonlat: upper right hand corner [longitude, latitude]
    """
    # convert to UTM60
    LLC_x,LLC_y = lonlat_utm(lon_or_x=LLC_lonlat[0],lat_or_y=LLC_lonlat[1])
    URC_x,URC_y = lonlat_utm(lon_or_x=URC_lonlat[0],lat_or_y=URC_lonlat[1])
    
    print("+LLHC\n\tLON/LAT: ({lon} ,{lat})\n\tUTM-60 X/Y: ({x} ,{y})".format(
                                            lon=LLC_lonlat[0],lat=LLC_lonlat[1],
                                            x=LLC_x,y=LLC_y))
    print("+URHC\n\tLON/LAT: ({lon} ,{lat})\n\tUTM-60 X/Y: ({x} ,{y})".format(
                                            lon=URC_lonlat[0],lat=URC_lonlat[1],
                                            x=URC_x,y=URC_y))
    
    # flooring all UTM values to get even numbers
    LLC_x = int(LLC_x)
    LLC_y = int(LLC_y)
    URC_x = int(URC_x)
    URC_y = int(URC_y)

    # potentially add a section which changes the URC values to give clean
    # difference values in the next section
    
    # determine the current ratio given input lat lon values
    x_diff = abs(LLC_x-URC_x)
    y_diff = abs(LLC_y-URC_y)
    if x_diff > y_diff:
        xy = 'X/Y'
    else:
        xy = 'Y/X'

    raw_ratio = max([x_diff,y_diff])/min([x_diff,y_diff])
    
    # provide candidates for simple aspect ratios to compare against raw ratio
    difflist = []
    startindex, endindex = 2,10
    print("\nBetween {}/{}".format(startindex,startindex-1),end=" ")
    for numerator in range(startindex,endindex):
        denominator = numerator - 1
        comparison_ratio = numerator/denominator
        difference = abs(raw_ratio-comparison_ratio)
        difflist.append(difference)
    print("and {}/{}".format(numerator,denominator),end=", ")
    
    diffindex = difflist.index(min(difflist))
    int_numer = diffindex+startindex
    int_denom = diffindex+startindex-1
    new_ratio = int_numer/int_denom
    print('ratio is closest to {num}/{den}'.format(num=int_numer,
                                                            den=int_denom))
    
    # suggest different lat/lon values to get exact aspect ratio
    # preferentially move upper edge
    if xy == 'X/Y':
        new_ratio = 1/new_ratio
    
    URC_y_new = (URC_x - LLC_x)*new_ratio + LLC_y  
    print('Moving upper edge by {}m to compensate\n'.format(
                                                        abs(URC_y_new-URC_y)))
                                                        
    # print new coordinates and final grid spacings for given resolutions
    URC_lon_new,URC_lat_new = lonlat_utm(lon_or_x=URC_x,
                                         lat_or_y=URC_y_new,
                                         inverse=True)
                                         
    print("NEW COORDINATES (WGS-84):")
    outstr=("LLC (lon,lat):\t{0},{1}\n"
            "URC (lon,lat):\t{2:.3f},{3:.3f}\n")
    print(outstr.format(LLC_lonlat[0],LLC_lonlat[1],URC_lon_new,URC_lat_new))
                                             
    print("NEW COORDINATES (UTM-60):")
    outstr=("longitude_min:\t{:.1f}\nlongitude_max:\t{:.1f}\n"
            "latitude_min:\t{:.1f}\nlatitude_max:\t{:.1f}\n")
    print(outstr.format(LLC_x,URC_x,LLC_y,URC_y_new))
    

    
    dx_out = (URC_x - LLC_x) * 1E-3
    dy_out = (URC_y_new - LLC_y) * 1E-3
    
    print("the region is {:.1f} km by {:.1f} km\n".format(dx_out,dy_out))
    # given grid spacing s, suggest values for NX and NY
    for s in [0.25,0.5,0.75,1,4,8]:
        print("NX={}, NY={} for {}KM SPACING \n".format(
                                            int((URC_x-LLC_x)*1/s*1E-3),
                                            int((URC_y_new-LLC_y)*1/s*1E-3),s)
                                            )
                                            
    
# ================================== EXAMPLES ==================================
def example_call_lonlat_utm(choice=0):
    """example call for function lonlat_utm, different initial coordinate calls
    separated by commented sections
    """
    # original kaikoura to east cape
    if choice == 0:
        LLC_lonlat = (173.3795,-42.4198)
        URC_lonlat = (178.8848,-37.4556)
        LLC_xy= lonlat_utm(LLC_lonlat[0],LLC_lonlat[1])
        URC_xy= lonlat_utm(URC_lonlat[0],URC_lonlat[1])
        print(LLC_xy)
        print(URC_xy)
	
	#modified kaikoura to east cape
    elif choice == 1:
        LLC_lonlat = (171311,5286948)
        URC_lonlat = [LLC_lonlat[0]+495,LLC_lonlat[1]+566]
        LLC_xy= lonlat_utm(LLC_lonlat[0],LLC_lonlat[1],inverse=True)
        URC_xy= lonlat_utm(URC_lonlat[0],URC_lonlat[1],inverse=True)
        print(LLC_xy)
        print(URC_xy)

    # carls project
    elif choice == 2:
        TAPE_LLC = (44429.4,5335568.8)
        TAPE_URC = (714429.4,5835568.8)
        LLC_xy= lonlat_utm(TAPE_LLC[0],TAPE_LLC[1],inverse=True)
        URC_xy= lonlat_utm(TAPE_URC[0],TAPE_URC[1],inverse=True)
        print(LLC_xy)
        print(URC_xy)
    
def example_call_determine_nxny():
    """example call for determine_nxny
    """
    # kaikoura to east cape (changed from -42.5 LLC lat 28.8.18)
    LLC_lonlat = (170,-42.8)
    URC_lonlat = (179.5,-37)
    determine_nxny(LLC_lonlat,URC_lonlat)
    
if __name__ == "__main__":
    example_call_determine_nxny()


    
    
    
