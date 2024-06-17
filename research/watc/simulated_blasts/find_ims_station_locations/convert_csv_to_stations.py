"""
Convert the Google Sheets .csv file into a STATIONS file to get "theoretical"
station locations for synthetic testing
"""
import numpy as np

net = "IM"
elv = 0.0
bur = 0.0
fmt = "{sta:>6}{net:>6}{lat:12.4f}{lon:12.4f}{elv:7.1f}{bur:7.1f}\n"

primary = np.loadtxt("ims_primary_seismic.csv", skiprows=1, usecols=[0,1,2,3], 
                     delimiter=",", dtype="object")

auxiliary = np.loadtxt("ims_auxiliary_seismic.csv", skiprows=1, 
                       usecols=[0,1,2,3], delimiter=",", dtype="object")

fid = "/Users/chow/Repos/spectral/research/watc/simulated_blasts/DATA/STATIONS_IMS_PLANNED_LOCATIONS"
with open(fid, "w") as f:
    for p in primary:
        sta, lat, lon, _ = p
        f.write(fmt.format(sta=sta, net=net, lat=float(lat), lon=float(lon), 
                           elv=elv, bur=bur))
    for a in auxiliary:
        sta, lat, lon, _ = a
        f.write(fmt.format(sta=sta, net=net, lat=float(lat), lon=float(lon), 
                           elv=elv, bur=bur))
