"""script for looking at the hobitss data
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from getdata import pathnames, event_stream
from obspy.clients.fdsn import Client
from obspy import read, read_inventory, read_events

import matplotlib as mpl
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
print("===========================================")
print("==== Be aware: ignoring FutureWarnings ====")
print("===========================================")


def download_hobitss_data(event):
    """download hobitss data via fdsn using obspy event object
    """
    # determine start and end of hobitss experiment
    inv = read2014p864702_LOBS_inventory(
                                    './datafiles/hobitss/hobitss_stations.xml')
    starts,ends = [],[]
    for sta in inv[0]:
        starts.append(sta[0].start_date)
        ends.append(sta[0].end_date)
    global_start = max(starts)
    global_end = min(ends)

    # event information for waveform gather
    eventid = str(event.resource_id).split('/')[1]
    starttime = event.origins[0].time - 100
    endtime = starttime + 1500

    # skip any event that falls outside the experiment
    if (starttime <= global_start) or (endtime >= global_end):
        print(eventid,"not within timeframe")
        return

    print(event)
    c = Client("IRIS")
    EBS_data = c.get_waveforms(network="YH",
                               station="EBS*",
                               location="",
                               channel="*",
                               starttime=starttime,
                               endtime=endtime,
                               attach_response=True)
    LOBS_data = c.get_waveforms(network="YH",
                               station="LOBS*",
                               location="",
                               channel="HH*",
                               starttime=starttime,
                               endtime=endtime,
                               attach_response=True)

    # write waveform data
    EBS_data.write('./datafiles/hobitss/{}_EBS.mseed'.format(eventid),
                                                                format='MSEED')
    LOBS_data.write('./datafiles/hobitss/{}_LOBS.mseed'.format(eventid),
                                                                format='MSEED')

def preprocess(st,inv):
    """preprocess waveform data
    """
    st_manipulate = st.copy()
    st_manipulate.detrend("demean")
    st_manipulate.detrend("linear")
    st_manipulate.taper(max_percentage=0.05)
    st_manipulate.attach_response(inv)
    st_manipulate.remove_response(output="VEL",water_level=0)

    return st_manipulate

def pretty_grids(input_ax):
    """make dem grids pretty
    """
    input_ax.set_axisbelow(True)
    input_ax.tick_params(which='both',direction='in',top=True,right=True)
    input_ax.minorticks_on()
    input_ax.grid(which='minor',
                    linestyle=':',
                    linewidth='0.5',
                    color='k',
                    alpha=0.25)
    input_ax.grid(which='major',
                    linestyle='-',
                    linewidth='0.5',
                    color='k',
                    alpha=0.15)

def align_yaxis(ax1,v1,ax2,v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1
    """
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)

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
        ax.ticklabel_format(style='sci',
                                axis='y',
                                scilimits=(0,0))
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

def amplitude_threshold(t,tr,threshold_percentage):
    """based on some amplitude threshold, return all points that fall above,
    and the time length covered
    """
    duration_a,tover_plot,aover_plot,threshold_plot,sample_plot = [],[],[],[],[]

    samp_rate = tr.stats.sampling_rate
    tr_manipulate = np.sqrt(tr.data**2)
    peak_amp = tr_manipulate.max()
    threshold = peak_amp * threshold_percentage
    threshold_plot.append(threshold)

    # loop over seismogram, determine start and end of peak energy
    # a for amplitude, s for sample
    a_over, s_over = [],[]
    for i,amplitude in enumerate(tr_manipulate):
        if amplitude >= threshold:
            s_over.append(i)
            a_over.append(amplitude)

    # find edgepoints by checking if the next sample j is the same as i+1
    s_edge,a_edge,sections,samples = [s_over[0]],[a_over[0]],[],[]
    for i,(S,A) in enumerate(zip(s_over[1:-2],a_over[1:-2])):
        if s_over[i+2] != (S + 1):
            # determine number of samples covered
            samples.append(len(s_edge))
            s_edge,a_edge = [s_over[i+1]],[a_over[i+1]]
        else:
            # if the next sample is the same, keep going
            s_edge.append(S)
            a_edge.append(A)

    duration_in_seconds = round(sum(samples)/samp_rate,0)

    # convert samples to time
    t_over = []
    for S in s_over:
        t_over.append(t[S])

    return t_over,a_over,duration_in_seconds

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
# for i in range(1,10):
#     station_select = "LOBS{}".format(i)
#     print(station_select)
#     try:
#         lobs = all_LOBS.select(station=station_select)
#         st = preprocess(lobs,inv)
#         st.filter('bandpass',freqmin=1/30,freqmax=1/5)
#         f = component_subplots(st,lobs)
#         figurename = "{e}_{i}.png".format(e=eventid,i=station_select)
#         plt.savefig(figure_out_path + figurename)
#         # plt.show()
#     except Exception as e:
#         print(e)
#         continue

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
