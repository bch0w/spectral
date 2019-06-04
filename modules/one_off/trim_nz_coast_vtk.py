"""The vtk file that Yoshi gave me contains the entire country, trim
that down to fit the dimensions of my mesh
"""
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

map_corners = [-42.65, -36.8, 172.95, 179.5]
nx_min, ny_min = lonlat_utm(map_corners[2], map_corners[0])
nx_max, ny_max = lonlat_utm(map_corners[3], map_corners[1])

coast_vtk = "/Users/chowbr/Documents/subduction/data/VTK/nz_resf_coast_mod_utm60H_xyz.vtk"
with open(coast_vtk, "r") as f:
    lines = f.readlines()

lines_out = lines[:5]
print("original file: {}".format(len(lines)))
for i, line in enumerate(lines[5:]):
    if "LINES" in line:
        for line_final in lines[i:]:
            lines_out.append(line_final)
        break
    nx, ny = line.split()[:2]
    if (nx_min < float(nx) < nx_max) and (ny_min < float(ny) < ny_max):
        lines_out.append(line)

print("cut file: {}".format(len(lines_out)))
coast_vtk_out = "/Users/chowbr/Documents/subduction/data/VTK/nz_resf_coast_mod_utm60H_xyz_trimmed.vtk"
import ipdb;ipdb.set_trace()
with open(coast_vtk_out, "w") as f:
    for line in lines_out:
        f.write(line)


        
    

