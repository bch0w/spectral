"""
Plot all MTs with PyGMT
"""
import pygmt
import numpy as np
from glob import glob

region = [0, 16, 0, 1]
projection = "X10c/4c"
frame = ["af", "+ggray90"]
fig = pygmt.Figure()
fig.basemap(region=region, projection=projection, frame=frame)

for i, fid in enumerate(sorted(glob("CMTSOLUTION_*"))):
    mrr, mtt, mpp, mrt, mrp, mtp = np.loadtxt(fid, skiprows=7, usecols=1, 
                                              dtype=float)
    mt = {"mrr": mrr, "mtt": mtt, "mff": mpp, 
          "mrt": mrt, "mrf": mrp, "mtf": mtp, 
          "exponent": 1}

    fig.meca(spec=mt, scale="0.5c", longitude=i+1, latitude=0.5, depth=0)
fig.savefig("lunes.png")


