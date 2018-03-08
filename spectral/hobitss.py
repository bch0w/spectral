"""script for looking at the hobitss data
"""
import sys
sys.path.append('../modules')
import os
import numpy as np
import matplotlib.pyplot as plt
from obspy import read, read_inventory, read_events

# module functions
from getdata import pathnames, event_stream
from plotmod import pretty_grids, align_yaxis
from procmod import preprocess, amplitude_threshold

import matplotlib as mpl
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
print("===========================================")
print("==== Be aware: ignoring FutureWarnings ====")
print("===========================================")


def component_subplots(st,st_raw=False):
    """subplots of 3 components for a stream of a single instrument,
    also the option of plotting raw seismogram in the background, normalized
    by the new data
    """
    # some initial preprocessing
    stats = st[0].stats
    t = np.linspace(-100,1500,stats.npts)

    # set plot
    f,(ax1,ax2,ax3) = plt.subplots(3,figsize=(9,5),
                                    sharex=True,
                                    sharey=False,
                                    dpi=200)
    # plot st
    axes = [ax1,ax2,ax3]
    for tr,ax in zip(st,axes):
        ax.plot(t,tr.data,'k',label=tr.get_id(),zorder=3)
        tover,aover,duration = amplitude_threshold(t,tr,0.15)
        ax.scatter(tover,np.zeros(len(tover)),
                        c='r',marker='x',s=0.05,zorder=4,
                        label="duration {} s".format(duration))
        ax.legend(prop={'size':5})
        pretty_grids(ax)

    # plot raw if available
    if st_raw:
        ax1a,ax2a,ax3a = ax1.twinx(),ax2.twinx(),ax3.twinx()
        twin_axes = [ax1a,ax2a,ax3a]
        for tr,tax in zip(st_raw,twin_axes):
            tax.plot(t,tr.data,'gray',alpha=0.25,zorder=2)
            tax.ticklabel_format(style='sci',
                                    axis='y',
                                    scilimits=(0,0))

    # plot manipulation
    ax1.set_title("2014p864702")
    ax2.set_ylabel("velocity (m/s)")
    ax2a.set_ylabel("raw (counts)")
    ax3.set_xlabel('time (sec)')
    plt.xlim([-50,1000])

    return f

# =================================== MAIN ====================================
# eventid = "2014p715167"
eventid = "2014p864702"
datafiles = pathnames()["hobitss"]
cat = read_events(os.path.join(datafiles,"hobitss_events.xml"))
inv = read_inventory(os.path.join(datafiles,"hobitss_stations.xml"))
all_LOBS = read(os.path.join(datafiles,"{}_LOBS.mseed".format(eventid)))
all_EBS = read(os.path.join(datafiles,"{}_EBS.mseed".format(eventid)))

# output files
plotspath = pathnames()["plots"]
figure_out_path = os.path.join(plotspath,"hobitss",'')

# LOBS DATA
for i in range(1,10):
    station_select = "LOBS{}".format(i)
    print(station_select)
    try:
        lobs = all_LOBS.select(station=station_select)
        st = preprocess(lobs,inv)
        st.filter('bandpass',freqmin=1/30,freqmax=1/5)
        f = component_subplots(st,lobs)
        figurename = "{e}_{i}.png".format(e=eventid,i=station_select)
        # plt.savefig(figure_out_path + figurename)
        plt.show()
    except Exception as e:
        print(e)
        continue

# ONLAND DATA
# station_list = ["KNZ","BFZ","BKZ","TSZ","PUZ"]
station_list = ['WAZ']
for station_choice in station_list:
    st_gn,inv_gn,_ = event_stream(station=station_choice,
                                channel="HH*",
                                event_id=eventid,
                                startpad=100,
                                endpad=1000)
    st = preprocess(st_gn,inv_gn)
    st.filter('bandpass',freqmin=1/30,freqmax=1/5)
    f = component_subplots(st,st_gn)
    figurename = "{e}_{i}.png".format(e=eventid,i=station_choice)


    plt.savefig(figure_out_path + figurename)
