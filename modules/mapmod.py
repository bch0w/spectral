"""module primarily for generating basemap objects of New Zealand for plotting,
for example, station locations, event locations, event beachballs etc.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm

import matplotlib as mpl
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1

from getdata import pathnames
from procmod import myround

def inventory_to_npz_array(filepath):
    """for quick plotting, generate npz arrays containing latitude,longitude
    and names of stations
    """
    from obspy import read_inventory
    inv = read_inventory(filepath)
    network = inv[0].code
    lats,lons,names = [],[],[]
    for sta_ in inv[0]:
        station_name = "{N}.{S}".format(N=network,S=sta_.code)
        names.append(station_name)
        lats.append(sta_.latitude)
        lons.append(sta_.longitude)

    new_filename = os.path.splitext(filepath)[0]
    np.savez(new_filename,NAME=names,LAT=lats,LON=lons)

def pd_pickle_to_npz_array(filepath):
    import pandas as pd
    cat = pd.read_pickle(filepath)
    names,lats,lons,mags = [],[],[],[]
    for i_ in range(len(cat)):
        event = cat.iloc[i_]
        names.append(event.event_id)
        lats.append(event.latitude)
        lons.append(event.longitude)
        mags.append(event.mw)

    new_filename = os.path.splitext(filepath)[0] + '.npz'
    np.savez(new_filename,NAME=names,LAT=lats,LON=lons,MAG=mags)


def initiate_basemap(etopo_else_flat_map=False,map_corners=[-50,-32.5,165,180],
                        draw_lines=True,figsize=(10,10),dpi=100):
    """set up local map of NZ to be filled
    default map corners give a rough box around new zealand
    etopo_else_flat_map: bool to control continent fill with either low res topo
    or solid line representation
    """
    # flat map configuration
    continent_color = 'w'
    lake_color = 'w'

    # initiate map
    fig = plt.figure(figsize=figsize,dpi=dpi)
    m = Basemap(projection = 'stere',
                resolution = 'h',
                rsphere = 6371200,
                lat_0 = np.mean(map_corners[:2]),
                lon_0 = np.mean(map_corners[2:]),
                llcrnrlat = map_corners[0],
                llcrnrlon = map_corners[2],
                urcrnrlat = map_corners[1],
                urcrnrlon = map_corners[3],
                )
    if etopo_else_flat_map:
        m.etopo()
    else:
        m.fillcontinents(color=continent_color,lake_color=lake_color)
    m.drawcoastlines(linewidth=1.25)

    # draw parallels and meridians
    if draw_lines:
        step_by = 2.5
        start_parallels = myround(map_corners[0]) - step_by
        end_parallels = myround(map_corners[1]) + step_by
        start_meridians = myround(map_corners[2]) - step_by
        end_meridians = myround(map_corners[3]) + step_by

        parallels = np.arange(start_parallels,end_parallels,step_by)
        meridians = np.arange(start_meridians,end_meridians,step_by)

        m.drawparallels(parallels,labels=[1,1,1,1],linewidth=0)
        m.drawmeridians(meridians,labels=[1,1,1,1],linewidth=0)

    return fig, m

def populate_basemap(m,lats,lons,names=None,mags=None,cs='b',ce='orange'):
    """fill map with latitude/longitude pairs, i.e. stations, events
    """
    if names is None: names = []
    if mags is None: mags = []

    X,Y = m(lons,lats)
    # plotting stations
    if len(mags) == 0:
        scatter = m.scatter(X,Y,
                            marker='v',
                            color=cs,
                            s=50,
                            edgecolor='k',
                            zorder=5)
    # plotting earthquakes, scaled by size
    else:
        scatter = m.scatter(X,Y,
                            marker='o',
                            color=ce,
                            edgecolor='k',
                            s=mags**2,
                            zorder=6)
        mags = np.array([str(_) for _ in mags])
        for m_,x_,y_ in zip(mags,X,Y):
            plt.annotate(m_,xy=(x_,y_),xytext=(x_,y_),fontsize=5)

    if len(names) != 0:
        for n_,x_,y_ in zip(names,X,Y):
            plt.annotate(n_,xy=(x_,y_),xytext=(x_,y_),zorder=6,fontsize=8.5)


def generate_map():
    """full function to initiate and populate basemap
    """
    allstation = False
    map_corner_dict = {"NZ":[-50,-32.5,165,180],
                     "NORTHISLAND":[-42,-34,172,179],
                     "SOUTHISLAND":[-47.5,-40,165,175]}

    f,m = initiate_basemap(map_corners=map_corner_dict["NZ"])


    # plot stations
    stationfile = 'NZ_NI_BB_seismometers.npz'
    if allstation:
        stationfile = 'NZ_NI_ALL_seismometers.npz'
    npzfolder = pathnames()['data'] + 'STATIONXML'
    sta = np.load(os.path.join(npzfolder,stationfile))
    names,lats,lons = sta['NAME'],sta['LAT'],sta['LON']
    populate_basemap(m,lats,lons,names)

    # plot events
    eventfile = 'tomCat.npz'
    npzfolder = pathnames()['data'] + 'CATBUILD'
    evt = np.load(os.path.join(npzfolder,eventfile))
    names,lats,lons,mags = evt['NAME'],evt['LAT'],evt['LON'],evt['MAG']

    # filter events by magnitude
    indices = np.where(mags>5.0)[0]
    names = names[indices]
    lats = lats[indices]
    lons = lons[indices]
    mags = mags[indices]
    print(len(mags))
    # populate_basemap(m,lats,lons,names,mags)

    # filter events by event name
    names,lats,lons,mags = evt['NAME'],evt['LAT'],evt['LON'],evt['MAG']

    shortlist = ['2403682','2593170','2015p822263m','2015p768477','2799448',
                '2354133','2016p860224','2016p669820','2013p614135',
                '2014p240655','2016p881118','2016p860224']
    indices = np.array([],dtype='int64')
    for ev_ in shortlist:
        index = np.where(names==ev_)[0]
        indices = np.concatenate((indices,index),axis=0)
    names = names[indices]
    lats = lats[indices]
    lons = lons[indices]
    mags = mags[indices]
    populate_basemap(m,lats,lons,names,mags,ce='green')


    plt.show()


if __name__ == "__main__":
    generate_map()
