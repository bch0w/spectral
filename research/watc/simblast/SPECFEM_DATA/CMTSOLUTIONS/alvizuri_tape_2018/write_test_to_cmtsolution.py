"""
Write a CMTSOLUTION file from Celso's 2018 SRL paper
"""
from glob import glob
import numpy as np



events = np.loadtxt("alvizuri_tape_2018_supp_tab2.txt", skiprows=3, dtype=float)
for i, event in enumerate(events):
    y, m, d, h, mi, sec, lon, lat, depth, \
            mrr, mtt, mpp, mrt, mrp, mtp, m0, mw = event
    cnvt = 1E7  # convert N*m -> dyne*cm
    name = f"NK{i+1}"
    out = f"""\
 PDE {y.astype(int)}  {m.astype(int)}  {d.astype(int)}  {h.astype(int)} {mi.astype(int)} {sec:.2f} {lat} {lon} {depth} {mw} {mw} {name}
event name:     {name}
time shift:     00.0000
half duration:  00.0000
latitude:  {lat:12.4f}
longitude: {lon:12.4f}
depth:     {depth:12.4f}
Mrr:      {mrr * cnvt:13.6e} 
Mtt:      {mtt * cnvt:13.6e} 
Mpp:      {mpp * cnvt:13.6e} 
Mrt:      {mrt * cnvt:13.6e} 
Mrp:      {mrp * cnvt:13.6e} 
Mtp:      {mtp * cnvt:13.6e}"""

    with open(f"CMTSOLUTION_ALVIZURI2018_{name}", "w") as f:
        f.write(out)

