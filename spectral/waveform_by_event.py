"""Plot waveforms for a given event
"""
import os
import sys
import glob
import traceback
sys.path.append('../modules/')
import numpy as np
from os.path import join
from obspy import UTCDateTime, read, Stream
from obspy.geodetics.base import gps2dist_azimuth

# module functions
import getdata
import procmod
import plotmod
from getdata import pathnames

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=mpl.cbook.mplDeprecation)

# =================================== FUNC ====================================
def plot_streams(st,bounds,save=False,show=True):
    """plot 6 streams in 3 subplot figure
    """
    c1,c2,c3 = "Z","N","E"
    f,(ax1,ax2,ax3) = plt.subplots(3,sharex=True,sharey=False,
                                                    figsize=(9,5),dpi=200)
    # create time axis
    stats = st[0].stats
    net,sta,loc,cha = st[0].get_id().split('.')
    t = np.linspace(0,stats.endtime-stats.starttime,stats.npts)

    # axes lists for plotting
    axes = [ax1,ax2,ax3]
    obs_list = ["HH{}".format(c1),"HH{}".format(c2),"HH{}".format(c3)]

    # plot
    for ax,obs in zip(axes,obs_list):
        obs_select = st.select(channel=obs)[0]
        A = ax.plot(t,obs_select.data,color='k',label=obs_select.get_id())
        plotmod.pretty_grids(ax)

    # final plot
    ax1.set_xlim([t.min(),t.max()])
    ax1.set_title("{e} {s} [{t0}-{t1}s]".format(e=event_id,
                                                    s=stats.station,
                                                    t0=bounds[0],
                                                    t1=bounds[1])
                                                    )

    ax1.set_ylabel(c1)
    ax2.set_ylabel("velocity (m/s)\n{}".format(c2))
    ax3.set_ylabel(c3)
    ax3.set_xlabel("time (sec)")

    if save:
        figtitle = "{e}_{s}.png".format(e=event_id,s=stats.station)
        figfolder = join(pathnames()['kupeplots'],"obsynth_plots",figtitle)
        plt.savefig(figfolder,dpi=200)

    if show:
        plt.show()

    return f

def plot_two_streams(st,bounds,twinax=False,save=False,show=True):
    """plot 6 streams in 3 subplot figure
    """
    c1,c2,c3 = "Z","N","E"
    f,(ax1,ax2,ax3) = plt.subplots(3,sharex=True,sharey=False,
                                                    figsize=(9,5),dpi=200)
    # create time axis
    stats = st[0].stats
    net,sta,loc,cha = st[0].get_id().split('.')
    t = np.linspace(0,stats.endtime-stats.starttime,stats.npts)

    # axes lists for plotting
    axes = [ax1,ax2,ax3]
    twin_axes = [ax1,ax2,ax3]
    if twinax:
        ax1a,ax2a,ax3a = ax1.twinx(),ax2.twinx(),ax3.twinx()
        twin_axes = [ax1a,ax2a,ax3a]
        ax2a.set_ylabel("velocity (m/s)")

    obs_list = ["HH{}".format(c1),"HH{}".format(c2),"HH{}".format(c3)]
    st_1 = st.select(station="*Z")
    st_2 = st.select(station="GA*")

    # plot
    for ax,tax,obs in zip(axes,twin_axes,obs_list):
        obs_select = st_1.select(channel=obs)[0]
        obs_2_select = st_2.select(channel=obs)[0]
        A = ax.plot(t,obs_select.data,color='k',label=obs_select.get_id())
        B = tax.plot(t,obs_2_select.data,color='r',label=obs_2_select.get_id())
        plotmod.pretty_grids(ax)
        lines = A+B
        labels = [l.get_label() for l in lines]
        ax.legend(lines,labels,prop={"size":5})
        plotmod.pretty_grids(tax)
        plotmod.align_yaxis(ax,0,tax,0)

    # final plot
    ax1.set_xlim([t.min(),t.max()])
    ax1.set_title("{e} {s} [{t0}-{t1}s]".format(e=event_id,
                                                    s=stats.station,
                                                    t0=bounds[0],
                                                    t1=bounds[1])
                                                    )

    ax1.set_ylabel(c1)
    ax2.set_ylabel("displacement (m)\n{}".format(c2))
    ax3.set_ylabel(c3)

    ax3.set_xlabel("time (sec)")

    if save:
        figtitle = "{e}_{s}.png".format(e=event_id,s=stats.station)
        figfolder = join(pathnames()['kupeplots'],"obsynth_plots",figtitle)
        plt.savefig(figfolder,dpi=200)

    if show:
        plt.show()


    return f

def initial_data_gather(code,event_id,bounds):
    """grab data
    """
    st,inv,cat = getdata.event_stream(code=code,
                                        event_id=event_id)
    # bannister data
    st_ban = Stream()
    for Co in ["HHZ","HHN","HHE"]:
        ban_path = glob.glob(pathnames()['mseed'] + 'GA/{e}/*GA01.{c}*'.format(
                                                            e=event_id,
                                                            c=Co))
        st_ban += read(ban_path[0])

    origin_time = cat[0].origins[0].time
    st_proc = procmod.preprocess(st,inv=inv,output="DISP")
    for tr in st_proc:
        tr.data *= (10**9) # meters to micrometers
    st_ban_proc = procmod.preprocess(st_ban)

    st_alltogether = st_proc + st_ban_proc
    # st_alltogether.trim(origin_time,origin_time+300)
    procmod.trimstreams(st_alltogether)
    st_alltogether.filter("bandpass",freqmin=1/bounds[1],freqmax=1/bounds[0])


    return st_alltogether

# =================================== MAIN ====================================
event_list = glob.glob(pathnames()['mseed'] + 'GA/*')
code = "NZ.KNZ.*.HH*"
# event_id = "2013p545288"
bounds = [6,30]

check_list = []
event_work = []
for event in event_list:
    try:
        event_id = os.path.basename(event)
        st = initial_data_gather(code,event_id,bounds)
        plot_two_streams(st,bounds=bounds,twinax=False,show=True,save=False)
    except Exception as e:
        print("="*79)
        traceback.print_exc()
        print("="*79)
        continue
    check = input("y or n")
    check_list.append(check)
    event_work.append(event)

for i,j in zip(event_work,check_list):
    print(os.path.basename(i),j)
import ipdb;ipdb.set_trace()
