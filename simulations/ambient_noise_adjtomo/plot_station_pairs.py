"""
Plot the station pairs with source receiver connecter lines to see 
covereage over Alaska

https://stackoverflow.com/questions/21352580/matplotlib-plotting-numerous-disconnected-line-segments-with-different-colors
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections as mc


# stations = np.loadtxt("STATIONS_FILTERED", dtype=str)
stations = np.loadtxt("STATIONS_LIU2022", dtype=str)
stadict = {}
for station in stations:
    stadict[f"{station[1]}.{station[0]}"] = (float(station[2]), 
                                             float(station[3]))

with open("station_pairs.txt", "r") as f:
    lines = f.readlines()
    print(len(lines))

# Takes a while but compiles all available station pairs to plot a line for
# Do this one time and write to a text file because it is too many points to
# store in ram 

station_pairs_fid = "plot_station_pairs.txt"
if not os.path.exists(station_pairs_fid):
    connected = []
    sta_a = None
    with open(station_pairs_fid, "w") as f:
        for i, line in enumerate(lines):
            if len(line.split()) > 1:
                _, _, sta_a = line.strip().split()

                if sta_a not in stadict:
                    sta_a = None
                    continue

                lat_a, lon_a = stadict[sta_a]
                lon_a = lon_a % -360

            # Master stations not in the STATION list get skipped completely until
            # another master station comes up in the file
            if sta_a is None:
                continue

            sta_b = line.strip()
            if sta_b not in stadict:
                continue

            if f"{sta_a} {sta_b}" in connected or f"{sta_b} {sta_a}" in connected:
                continue

            lat_b, lon_b = stadict[sta_b]
            lon_b = lon_b % -360

            f.write(f"{lon_a} {lat_a} {lon_b} {lat_b}\n")
            connected.append(f"{sta_a} {sta_b}")

# Load in the coordinates
plot_lines = np.loadtxt(station_pairs_fid, dtype=float)
plot_tups = [[(a, b), (c, d)] for a, b, c, d in plot_lines]

# Get station markers to plot
lats = stations[:, 2].astype(float)
lons = stations[:, 3].astype(float)
lons = lons % -360
plt.scatter(lons, lats, c="None", edgecolor="g", linewidth=0.75, marker="v", 
            zorder=5)

# Plot connecting lines for all stations
lc = mc.LineCollection(plot_tups, colors="k", linewidth=0.5, alpha=0.01,
                       zorder=3)

ax = plt.gca()
ax.add_collection(lc)

ax.set_xlabel("Longitude (deg)")
ax.set_ylabel("Latitude (deg)")
ax.set_title(f"Liu et al. (2022) EGFs for Alaska;\n"
             f"STATIONS={len(stations)}; PAIRS={len(plot_lines)}")

plt.savefig("egf.png")
plt.show()

        
