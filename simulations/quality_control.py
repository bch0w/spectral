"""quality control to determine which events to use on tomography
"""
import sys
import traceback
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.append("../modules")
from obspy import read_inventory
from obspy.geodetics import locations2degrees

import getdata
import procmod
from getdata import pathnames

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


# =================================== PLOT =====================================
def mapper(tomCat):
    """26.4.18 outdated, use mapmod in modules from now on
    builds basemap using obspy and populates with list of lat/lon values
    """
    tomCat_by_mw = tomCat.sort_values('mw',ascending=False)
    lats,lons,infos = [],[],[]
    for i in range(50):
        event = tomCat_by_mw.iloc[i]
        event_id = event.event_id
        event_lat = event.latitude
        event_lon = event.longitude
        event_mw = event.mw
        event_info = "{}/{}".format(event_id,event_mw)
        lats.append(event_lat)
        lons.append(event_lon)
        infos.append(event_info)

    stationxml_file = (pathnames()['data'] +
            'STATIONXML/plotNZ_ALL_BB_seismometers_stations_only.xml')
    stations = read_inventory(stationxml_file)
    fig = stations.plot(marker='.',
                projection="local",
                resolution = "i",
                # color=color_dict[choice],
                size=25,
                continent_fill_color="white",
                water_fill_color="white",
                color_per_network=True,
                label=False,
                show=False)
    map_x,map_y = fig.bmap(lons,lats)
    scatter = fig.bmap.scatter(map_x,map_y,
                        marker='v',
                        s=150,
                        zorder=1000,
                        # c=colors,
                        edgecolor='k')
    for X,Y,M in zip(map_x,map_y,infos):
        plt.annotate(M,xy=(X,Y),xytext=(X,Y),
                fontsize=10,
                weight='bold')

    plt.show()
    sys.exit()

# =================================== FUNC =====================================
def nearby_stations(event_lat,event_lon,number_of=5):
    """find number_of closest geonet stations to a given lat/lon
    """
    coordinate_file = (pathnames()['data'] +
                    'STATIONXML/plotNZ_ALL_BB_seismometers_stations_only.npz')
    stationlist = np.load(coordinate_file)
    sta_lat = stationlist['LAT']
    sta_lon = stationlist['LON']
    sta_name = stationlist['NAME']
    distances = []
    for SLAT,SLON in zip(sta_lat,sta_lon):
        dist = locations2degrees(SLAT,SLON,event_lat,event_lon)[0]
        distances.append(dist)

    # sort by distances
    distlist = sorted(zip(sta_name,distances), key=lambda x:x[1])
    sta_name_sort, distances_sort = zip(*distlist)

    return sta_name_sort[:number_of], distances_sort[:number_of]

def process_events(event,plot=False):
    """given an event column from pd dataframe, process station information
    """
    stations,distances = nearby_stations(event.latitude,
                                         event.longitude,
                                         number_of=10)

    # plot event on each of the nearby stations
    code_template = "{s}.10.HH?"
    event_id = event.event_id.values[0]
    print(event_id)
    SNR_list = []
    i = 0
    for sta in stations:
        if i == 3:
            break
        try:
            code = code_template.format(s=sta)
            st,inv,cat = getdata.event_stream(code,event_id,startpad=100)
            st_proc = procmod.preprocess(st,inv=inv)
            st_proc.filter('bandpass',freqmin=1/30,freqmax=1/6)
            if plot:
                st_proc.plot()
            # SNR = procmod.signal_to_noise(st_proc[0].data)
            # print(sta,SNR)
            # SNR_list.append(SNR)
            i+=1
        except Exception as E:
            traceback.print_exc()
            continue

    return SNR_list

def event_picker():
    """sort out the pickle file
    """
    pickle_file = pathnames()['kupe'] + 'tomCat/tomCat'
    tomCat = pd.read_pickle(pickle_file)
    # mapper(tomCat)
    shortlist = ['2922302','2016p661375','2403682','2593170','2015p822263',
                 '2015p768477','2799448','2354133','2016p858260','2016p860224']
    shortlist = ['2016p669820','2013p614135','2014p240655','2016p881118']
    for event_id in shortlist:
        event = tomCat.loc[tomCat['event_id']==event_id]
        # try:
        process_events(event,plot=True)
        # except KeyboardInterrupt:
        #     sys.exit()
        # except Exception as e:
            # print("error")


def errorprint():
    print("="*79)
    traceback.print_exc()
    print("="*79)

if __name__ == "__main__":
    event_picker()
