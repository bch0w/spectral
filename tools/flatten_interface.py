"""
Interpolate plate interface values from SPECFEM formatted tomography files
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree

choice = "vs"

val_list = ["x", "y", "z", "vp", "vs", "rho", "qp", "qs"]
interface = np.loadtxt("williams_hikurangi_interface_utm60_nonan.xyz")

for tag in ["mantle", "crust", "shallow"]:
    fid = f"tomography_model_{tag}.xyz"
    print(fid)
    values = np.loadtxt(fid, dtype=float, skiprows=4)
    xyz = values[:, :3]

    # Query nearest neighbor only for points that fall in the applicable depths
    zmin = xyz[:, 2].min()
    zmax = xyz[:, 2].max()
    idx = np.where((interface[:, 2] > zmin) & (interface[:, 2] < zmax))
    inter_xyz = interface[idx]

    # Use SciPy KDTree to find nearest neighbor of tomo grid given the interface
    kdtree = KDTree(xyz)
    dist, point = kdtree.query(inter_xyz)
    print(f"{len(point)} points found")
    print(f"min dist: {dist.min()}")
    print(f"max dist: {dist.max()}")

    # Grab the tomo model value and plot on a 2D plane
    v_idx = val_list.index(choice)
    xyv = values[:, [0, 1, v_idx]][point]  # Vs
    np.save(tag, xyv)



