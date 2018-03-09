import os
import sys
sys.path.append('../modules/')
import numpy as np
from os.path import join
from obspy import UTCDateTime, read, Stream, read_events

# module functinos
from getdata import pathnames, event_stream, get_moment_tensor
from synmod import stf_convolve, get_GCMT_solution
from procmod import preprocess, trimstreams
from plotmod import pretty_grids, align_yaxis

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
print("===========================================")
print("==== Be aware: ignoring FutureWarnings ====")
print("===========================================")

# =================================== FUNC ====================================
def initial_data_gather(event_id,code,tmin,tmax):
    """gather event information, observation and synthetic traces,
    preprocess all traces accordingly and return one stream object with 6 traces
    """
    # station information
    network,station,location,channel = code.split('.')
    component = channel[-1]

    # event information
    MT = get_GCMT_solution(event_id)
    moment_tensor = MT.focal_mechanisms[0].moment_tensor
    half_duration = (moment_tensor.source_time_function['duration'])/2
    mt_geonet = get_moment_tensor(event_id)

    # event time from CMTSOLUTION used in simulation
    CMTpath = join(pathnames()['data'],
                    'KUPEDATA','CMTSOLUTIONS','{}CMTSOLUTION'.format(event_id))
    CMTSOL = read_events(CMTpath)
    syn_start = CMTSOL[0].origins[0].time

    # grab synthetic data locally
    syntheticdata_path = join(pathnames()['syns'],event_id,'')
    syntheticdata = Stream()
    for c in ["N","E","Z"]:
        syntheticdata_filename = "{n}.{s}.BX{co}.semv.mseed".format(n=network,
                                                                s=station,
                                                                co=c)
        syntheticdata += read(join(syntheticdata_path,syntheticdata_filename))

    # grab observation data
    observationdata,inv,cat = event_stream(station=station,
                                            channel=channel,
                                            event_id=event_id,
                                            startpad=0,
                                            endpad=350)

    # preprocessing, instrument response, STF convolution (synthetics)
    observationdata_proc = preprocess(observationdata,inv=inv,output='DISP')
    syntheticdata_preproc = preprocess(syntheticdata)
    syntheticdata_proc = stf_convolve(syntheticdata_preproc,half_duration)
    for tr in syntheticdata_proc:
        tr.stats.starttime = syn_start

    syntheticdata_proc.integrate()
    # combine, common sampling rate, filter, trim common time
    st_IDG = observationdata_proc + syntheticdata_proc
    st_IDG.resample(50)
    st_IDG.filter('bandpass',freqmin=1/tmax,freqmax=1/tmin)
    trimstreams(st_IDG)

    return st_IDG

def plot_obsynth(st):
    """plot streams in 3 subplot figure
    """
    f,(ax1,ax2,ax3) = plt.subplots(3,sharex=True,sharey=False,
                                                    figsize=(9,5),dpi=200)
    # create time axis
    stats = st[0].stats
    t = np.linspace(0,stats.endtime-stats.starttime,stats.npts)

    # axes lists for plotting
    ax1a,ax2a,ax3a = ax1.twinx(),ax2.twinx(),ax3.twinx()
    axes = [ax1,ax2,ax3]
    twin_axes = [ax1a,ax2a,ax3a]
    obs_list = ["HHN","HHE","HHZ"]
    syn_list = ["BXN","BXE","BXZ"]

    # plot
    for ax,tax,obs,syn in zip(axes,twin_axes,obs_list,syn_list):
        obs_select = st.select(channel=obs)[0]
        syn_select = st.select(channel=syn)[0]
        A = ax.plot(t,obs_select.data,color='k',label=obs_select.get_id())
        B = tax.plot(t,syn_select.data,color='r',label=syn_select.get_id())
        lines = A+B
        labels = [l.get_label() for l in lines]
        ax.legend(lines,labels,prop={"size":5})
        pretty_grids(ax)
        pretty_grids(tax)
        align_yaxis(ax,0,tax,0)

    # final plot
    ax1.set_title(event_id)
    ax2.set_ylabel('HH* velocity (m/s)')
    ax2a.set_ylabel('BX* velocity (m/s)')
    ax3.set_xlabel('time (sec)')
    plt.show()

    return f


# =================================== MAIN ====================================
station_code = sys.argv[1].upper()
event_id = "2014p240655"
code = "NZ.{}..HHZ".format(station_code)

st = initial_data_gather(event_id,code,tmin=5,tmax=30)
# st.plot(equal_scale=False)
fig = plot_obsynth(st)


# import ipdb;ipdb.set_trace()
