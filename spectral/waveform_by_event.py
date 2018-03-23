"""Plot waveforms for a given event
"""
import os
import sys
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

# =================================== MAIN ====================================
for i in [6,9,10]:
    code = "XX.RD{:0>2}.10.*".format(i)
    # event_id = "2018p130600"
    event_id = "2018p093947"
    bounds = [2,30]
    st,inv,cat = getdata.event_stream(code=code,
                                        event_id=event_id)
    # plotmod.plot_event_station(inv,cat)
    origin_time = cat[0].origins[0].time
    st.trim(origin_time,origin_time+300)
    st_proc = procmod.preprocess(st,inv=inv)
    st_proc.filter("bandpass",freqmin=1/bounds[1],freqmax=1/bounds[0])
    plot_streams(st_proc,bounds=bounds,show=True,save=False)
    # except Exception as e:
    #     continue
