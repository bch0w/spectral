"""
For visuazliation of coastline in Paraview
"""
import csv
import numpy as np


# Parameters
coastline_npy = "nz_resf_coast_mod_utm60H_xyz.npy"

# These are the bounds for Trelis 139_162 mesh
x_min = 1.73E5
x_max = 6.36E5
y_min = 5.27E6
y_max = 5.91E6

data = np.load(coastline_npy)

# Set the bounds if none are given
x_min = x_min or data[:, 0].min()
x_max = x_max or data[:, 0].max()
y_min = y_min or data[:, 1].min()
y_max = y_max or data[:, 1].max()

# Enforce the bounds
x_too_small = np.where(data[:, 0] < x_min)[0]
x_too_large = np.where(data[:, 0] > x_max)[0]
y_too_small = np.where(data[:, 1] < y_min)[0]
y_too_large = np.where(data[:, 1] > y_max)[0]

to_remove = np.unique(np.concatenate((x_too_small, x_too_large, y_too_small, 
                                      y_too_large), 0))

# delete in place
data = np.delete(data, to_remove, 0)

# Save as a csv file for easy read into Paraview
fid = "nz_coast_utm60_x{x:.0f}_y{y:.0f}.csv".format(x=(x_max-x_min)*1E-3,
                                                    y=(y_max-y_min)*1E-3
                                                    )

np.savetxt(fid, data, delimiter=",", fmt='%1.9e')

# Hacky way to stick a header onto the data
with open(fid, "r+") as f:
    content = f.read()
    f.seek(0, 0)
    f.write("X,Y,Z\n")
    f.write(content)

