"""
NIWA provides their bathymetry data in a single column of depth values. Need to
convert that to a grid and cut it to size

The NIWA grid is provided in EPSG: 3994
"""
import numpy as np
from pyproj import CRS, Transformer

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

# Define bounds of the NZNorth mesh in UTM60
xmin = 171311.85
xmax = 633468.31 
xmax += (xmax - xmin) / 3
ymin = 5286952.00
ymax = 5904085.00

# Generate the grid values based on the given header info
x_end = header["xllcorner"] + (header["ncols"] * header["cellsize"])
xvals = np.arange(header["xllcorner"], x_end, header["cellsize"])
assert(len(xvals) == header["ncols"])

y_end = header["yllcorner"] + (header["nrows"] * header["cellsize"])
yvals = np.arange(header["yllcorner"], y_end, header["cellsize"])
assert(len(yvals) == header["nrows"])

# Flatten the 2D grid into a long XYZ list
xgrid, ygrid = np.meshgrid(xvals, yvals)
x_vector = xgrid.ravel()
y_vector = ygrid.ravel()

# Convert the grid from EPGS:3994 to UTM60S
x_vector, y_vector = transform.transform(x_vector, y_vector)  # takes a while
data_vector = data_cut[::-1].ravel()
data_out = np.vstack((x_vector, y_vector, data_vector)).T

# Cut around the NZNorth Mesh bounds
idx = np.where((data_out[:, 0] > xmin) &
               (data_out[:, 0] < xmax) &
               (data_out[:, 1] > ymin) &
               (data_out[:, 1] < ymax))[0]
data_cut = data_out[idx]

np.savetxt("niwa_nznorthext_bathy_utm60s_250m.txt", data_cut)


