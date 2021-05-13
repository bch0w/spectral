"""
NIWA provides their bathymetry data in a single column of depth values. Need to
convert that to a grid and cut it to size

The NIWA grid is provided in EPSG: 3994
"""
import numpy as np
from pyproj import CRS, Transformer

def find_nearest(array, value):
    """
    Used to find nearest point in the grid based on value differences
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

with open("nzbathymetry_2016_ascii-grid.txt") as f:
    lines = f.readlines()

# Break the header into a dict object
header_lines = lines[:6]
header = {}
for line in header_lines:
    key, val = line.strip().split()
    header[key] = float(val)

try:
    data = np.load("nzbathymetry_2016_ascii-grid.npy")
except FileNotFoundError:
    data = lines[6:]
    for i, line in enumerate(data[:]):
        print(f"{i}/{len(data)}")
        # Convert to a numpy array
        data[i] = np.array(list(map(float, line.strip().split())))
    data = np.array(data)

# We will use PyProj to convert from EPSG:3994 to UTM60S (a.k.a EPSG:32760) 
crs_32760 = CRS.from_epsg(32760)  # UTM60S
crs_3994 = CRS.from_epsg(3994)
transform = Transformer.from_crs(crs_3994, crs_32760)
transform_r = Transformer.from_crs(crs_32760, crs_3994)

# Define bounds of the NZNorth mesh in UTM60
xmin_utm = 171312
xmax_utm = 633468
ymin_utm = 5286950
ymax_utm = 5890961

# Convert bounds to EPGS 3994
xmin, ymin = transform_r.transform(xmin_utm, ymin_utm)
xmax, ymax = transform_r.transform(xmax_utm, ymax_utm)
# xmin = 6141868.51 
# xmax = 6604612.03
# ymin = -3935671.97
# ymax = -3335636.67

# Generate the grid values based on the given header info
x_end = header["xllcorner"] + (header["ncols"] * header["cellsize"])
xvals = np.arange(header["xllcorner"], x_end, header["cellsize"])
assert(len(xvals) == header["ncols"])

y_end = header["yllcorner"] + (header["nrows"] * header["cellsize"])
yvals = np.arange(header["yllcorner"], y_end, header["cellsize"])
assert(len(yvals) == header["nrows"])

# Figure out what index values correspond to, manually choosing max x
ixmin = find_nearest(xvals, xmin)
ixmax = 8000 # find_nearest(xvals, xmax)

# Need to 'flip' the y-axis because we're dealing with a -Z axis
iymax = int(header["nrows"] - find_nearest(yvals, ymin))
iymin = int(header["nrows"] - find_nearest(yvals, ymax))

# Cut the data to the appropriate range. 
data_cut = data[iymin:iymax][:,ixmin:ixmax]
xvals = xvals[ixmin:ixmax]
yvals = yvals[iymin:iymax]

# Flatten the 2D grid into a long XYZ list
xgrid, ygrid = np.meshgrid(xvals, yvals)
x_vector = xgrid.ravel()
y_vector = ygrid.ravel()
x_vector, y_vector = transform.transform(x_vector, y_vector)
data_vector = data_cut[::-1].ravel()

data_out = np.vstack((x_vector, y_vector, data_vector)).T
np.savetxt("niwa_nznorth_bathymetry_utm60s_250m.txt", data_out)




