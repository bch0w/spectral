"""
Trying to extract isocontours from Williams et al. 2013 revised Hikurangi plate
interface model in UTM 60S coordinate system for use in Paraview plotting
"""
import numpy as np
from scipy.io import netcdf_file
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from pyatoa.utils.srcrcv import lonlat_utm


def get_contour_verts(cn):
    """
    From StackOverflow, get contour lines as an array of x,y vertices
    https://stackoverflow.com/questions/18304722/
        python-find-contour-lines-from-matplotlib-pyplot-contour
    """
    contours = []
    # for each contour line
    for cc in cn.collections:
        paths = []
        # for each separate section of the contour line
        for pp in cc.get_paths():
            xy = []
            # for each segment of that section
            for vv in pp.iter_segments():
                xy.append(vv[0])
            paths.append(np.vstack(xy))
        contours.append(paths)

    return contours

# Define the bounds of the output grid
x_min = 171312
x_max = 633468
y_min = 5286950
y_max = 5904080
z_val = None
cl = [-100, -75, -50, -40, -30, -20, -15, -12, -9, -6, -3]
grid_spacing = 100
fid_out = "contours_depth_utm60s.vtk"
            
# Read the .grd file from Williams
with netcdf_file("grid_exclude_wgs84.grd") as data:
    x = np.copy(data.variables["x"][:])
    y = np.copy(data.variables["y"][:])
    z = np.copy(data.variables["z"][:])

# Convert grid-like data to flat arrays
x, y = np.meshgrid(x, y)
x = x.flatten()
y = y.flatten()
z = z.flatten()

# Convert to UTM60S
x, y = lonlat_utm(x, y, -60, False)

# Create arrays for interpolation
points = np.vstack((x, y)).T
values = z

# Remove NaN values
idx = np.where(~np.isnan(z))
points = points[idx]
values = values[idx]

# Set up the regular grid
x_range = np.arange(x_min, x_max, grid_spacing)
y_range = np.arange(y_min, y_max, grid_spacing)
x_grid, y_grid = np.meshgrid(x_range, y_range)

# Interpolate interface onto regular grid
interp = griddata(points=points, values=values, xi=(x_grid, y_grid))

# Generate contour plot and extract vertices
cs = plt.contour(x_grid, y_grid, interp, levels=cl)
contours = get_contour_verts(cs)

# Collect nested arrays into lists 
x, y, z = [], [], []
with open("contours_depth_utm60s.txt" ,'w') as f:
    for i, contour in enumerate(contours):
        for segment in contour:
            for vals in segment:
                x.append(vals[0])
                y.append(vals[1])
                z.append(cl[i] * 1E3)

# Write values to a vtk file
npts = len(x)
with open(fid_out, "w") as f:
    # Write header
    f.write("# vtk DataFile Version 2.0\n"
            "loop\n"
            "ASCII\n"
            "DATASET UNSTRUCTURED_GRID\n"
            f"POINTS {npoints} float\n")

    # Write data
    for x_, y_, z_ in zip(x, y, z):
        f.write(f"{x_:15.2f} {y_:15.2f} {z_:15.2f}\n")

    # Write cell data
    f.write(f"CELLS {npoints} {npoints*2}\n")
    for i in range(npoints):
        f.write(f"1 {i}\n")

    # Write cell types
    f.write(f"CELL_TYPES {npoints}\n")
    for i in range(npoints):
        f.write(f"1 \n")

    # Write Z Scalar values
    f.write(f"\nPOINT_DATA {npoints}\n"
            f"SCALARS Z_Value float 1\n"
            f"LOOKUP_TABLE default\n")
    for z_ in z:
        f.write(f"{z_}\n")
