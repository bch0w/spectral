"""
Plot the station pairs with source receiver connecter lines to see 
covereage over Alaska

https://stackoverflow.com/questions/21352580/matplotlib-plotting-numerous-disconnected-line-segments-with-different-colors
"""
import os
import cartopy.crs as ccrs
import cartopy.feature as cf
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections as mc


stations = np.loadtxt("STATIONS", dtype=str)
coords = {}
for station in stations:
    # lat, lon
    coords[f"{station[1]}.{station[0]}"] = (float(station[2]), 
                                            float(station[3]))

# Get actual data pairs
with open("station_pairs.txt", "r") as f:
    lines = f.readlines()
srcrcvdict = {}
for line in lines:
    if not line.startswith("\t"):
        src_net, src_sta = line.strip().split()[-1].split(".")
        src = f"{src_net}.{src_sta}"
        srcrcvdict[src] = []
    else:
        srcrcvdict[src].append(line.strip())

# Initiate map
f = plt.figure(figsize=(11, 8), dpi=150)
ax = plt.axes(projection=ccrs.LambertConformal(central_latitude=(71.5+64.)/2, 
                                               central_longitude=(-169.5-137.5)/2)
              )
ax.set_extent([-169.5, -137.5, 64., 71.5])
ax.coastlines(lw=2)
ax.add_feature(cf.BORDERS)
ax.spines["geo"].set_visible(False)

# Loop through all source stations, plot
sources = np.loadtxt("SOURCES", dtype=str)
sta_plotted = [f"{_[1]}.{_[0]}" for _ in sources]
n = 0
for source in sources:
    src = f"{source[1]}.{source[0]}"

    src_lat, src_lon = coords[src]
    if src == "AK.D20K":
        ax.scatter(src_lon, src_lat, marker="o", edgecolor="k", linewidth=2, 
                   c="cyan",
                   s=90, zorder=20, transform=ccrs.PlateCarree())
    else:
        ax.scatter(src_lon, src_lat, marker="o", edgecolor="k", linewidth=2, 
                   c="salmon",
                   s=60, zorder=20, transform=ccrs.PlateCarree())
    if False:
        ax.text(src_lon-.15, src_lat+.175, src, transform=ccrs.PlateCarree(),
                c="k", fontsize=10, zorder=20, ha="right")

    # Plot all the stations connected 
    stations = srcrcvdict[src]
    for sta in stations:
        try:
            sta_lat, sta_lon = coords[sta]
        except KeyError:
            continue

        if False:
            ax.text(sta_lon+.2, sta_lat+.05, sta, transform=ccrs.PlateCarree(),
                    c="r", fontsize=10, zorder=20)

        # Plot the station marker if it hasn't shown up yet
        if sta not in sta_plotted:
            ax.scatter(sta_lon, sta_lat, marker="v", c="gold", 
                       linewidth=2, edgecolor="k",
                       s=60, zorder=19, transform=ccrs.PlateCarree())
            sta_plotted.append(sta)
        
        # Plot connecting lines between source and receiver

        if src == "AK.D20K":
            ax.plot([src_lon, sta_lon], [src_lat, sta_lat],
                    c="cyan", linewidth=2., alpha=0.5, transform=ccrs.Geodetic(),
                    zorder=19
                     )
        ax.plot([src_lon, sta_lon], [src_lat, sta_lat], 
                 c="k", linewidth=8., alpha=0.01, transform=ccrs.Geodetic(),
                 zorder=18
                 )
        n += 1

# plt.title(f"NSRC={len(sources)}; NRCV={len(stations)}; NPAIR={n}")
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,
                  zorder=10, linewidth=.5, alpha=0.2, color="k")
gl.xlabel_style= {"size": 12}
gl.ylabel_style= {"size": 12}
gl.top_labels=False
gl.right_labels=False
plt.savefig("srcrcv_pairs.png")
plt.show()

        
