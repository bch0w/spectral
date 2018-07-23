"""PARED AND EDITED SCRIPT TO MAKE EQ DURATION MAPS FOR VERTICAL COMPONENT
FIGURES FOR SLOW SLIP WORKSHOP PRESENTATION

use obspy to create station maps
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from getdata import pathnames
from obspy import read_inventory, read_events
from matplotlib.patches import Polygon

# ignore warnings
import warnings
warnings.filterwarnings("ignore",category=mpl.cbook.mplDeprecation)

def plot_polygon(fig):
    # plot low velocity overlay with a polygon
    nw_lat,nw_lon = -39.8764, 178.5385
    ne_lat,ne_lon = -39.0921, 177.1270
    se_lat,se_lon = -37.5630, 178.5248
    sw_lat,sw_lon = -38.3320, 179.9157
    x1,y1 = fig.bmap(nw_lon,nw_lat)
    x2,y2 = fig.bmap(ne_lon,ne_lat)
    x3,y3 = fig.bmap(se_lon,se_lat)
    x4,y4 = fig.bmap(sw_lon,sw_lat)
    poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],
                    facecolor='red',
                    edgecolor='k',
                    linewidth=3,
                    alpha=0.1)
    plt.gca().add_patch(poly)

# command line choices
event_list = ["2014p240655","2015p822263"]#,"2016p892721","2017p059122"]

# dictionaries for dynamic choices
component="vertical"

# grab station information
inv = read_inventory(pathnames()["spectral"] +
                        "duration/forSSW/north_island_geonet_stations.xml")


# grab event information
event_lat,event_lon,event_mag = [],[],[]
for event_id in event_list:
    event = read_events(pathnames()["spectral"] +
                                    "duration/forSSW/{}".format(event_id))[0]
    event_lat.append(event.origins[0].latitude)
    event_lon.append(event.origins[0].longitude)
    event_mag.append(round(event.magnitudes[0].mag,2))

# ========= manual plot scatterpoints - really hacked together =========
all_durations,all_stations,all_lats,all_lons = [],[],[],[]
for event_id in event_list:
    text_file = pathnames(
                )['spectral']+'duration/forSSW/{}_amplitudes.txt'.format(event_id)
    durations,durstations,lats,lons = [],[],[],[]
    with open(text_file,'r') as f:
        for line in f:
            eachline = line.split(' ')
            # skip over marked stations
            if "#" in eachline[0]:
                continue
            durstations.append(eachline[0])
            durations.append(float(eachline[1]))

    # grab latitude and longitude information
    for sta in durstations:
        sta_info = inv.select(station=sta)[0][0]
        lats.append(sta_info.latitude)
        lons.append(sta_info.longitude)

    all_durations.append(durations)
    all_stations.append(durstations)
    all_lats.append(lats)
    all_lons.append(lons)

# set colors by the maximum of lists
flat_durations = [item for sublist in all_durations for item in sublist]
D = [(_/max(flat_durations))**2 for _ in flat_durations]
cmap = plt.get_cmap('plasma')
flat_colors = cmap(D)
split_point = int(len(flat_durations)/2)
all_colors = [flat_colors[:split_point],flat_colors[split_point:]]

# for each event, plot a new map with the same colorbar
for LAT,LON,STA,DUR,COL,EV_LAT,EV_LON,EV_MAG,EID in zip(
                                    all_lats,all_lons,all_stations,
                                    all_durations,all_colors,event_lat,
                                    event_lon,event_mag,event_list):

    # plot stations
    fig = inv.plot(label=False,
            marker='.',
            projection="local",
            resolution = "i",
            color='r',
            size=6,
            show=False)


    plot_polygon(fig)

    map_x,map_y = fig.bmap(LON,LAT)
    # plot onto basemap
    scatter = fig.bmap.scatter(map_x,map_y,
                                marker='v',
                                s=150,
                                zorder=1000,
                                c=COL,
                                edgecolor='k',
                                linewidth=0.5)
    for i,(sta,dur) in enumerate(zip(STA,DUR)):
        dur = int(dur)
        plt.annotate("{}({}s)".format(sta,dur),
                        xy=(map_x[i],map_y[i]),
                        xytext=(map_x[i]+5000,map_y[i]+5000),
                        fontsize=12,
                        weight='bold',
                        zorder=201)

    # plot event
    eventx,eventy = fig.bmap(EV_LON,EV_LAT)
    fig.bmap.scatter(eventx,eventy,
                        marker='o',
                        s=150,
                        c='c',
                        zorder=200,
                        edgecolor='k')
    plt.annotate("M{}".format(EV_MAG),
                    xy=(eventx,eventy),
                    xytext=(eventx-7500*1.5,eventy+7500),
                    fontsize=12,
                    zorder=200,
                    weight='bold')

    plt.title("{}".format(EID))

    LLC = fig.bmap(173,-41)
    URC = fig.bmap(180,-37)
    plt.xlim([LLC[0],URC[0]])
    plt.ylim([LLC[1],URC[1]])


    plt.show()

# solo colorbar
fig = plt.figure(figsize=(8, 3))
ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
norm = mpl.colors.PowerNorm(vmin=20,
                            vmax=170,
                            gamma=2)
cb1 = mpl.colorbar.ColorbarBase(ax1,
                                cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
cb1.set_label("Shaking duration (seconds)")
plt.show()
