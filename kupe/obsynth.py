import os
import sys
sys.path.append('../modules/')
import numpy as np
from os.path import join
from obspy import UTCDateTime, read, Stream, read_events

# module functions
import getdata
import synmod
import procmod
import plotmod
from getdata import pathnames

# from getdata import pathnames, event_stream, \
#                     get_moment_tensor, get_GCMT_solution
# from synmod import stf_convolve
# from procmod import preprocess, trimstreams
# from plotmod import pretty_grids, align_yaxis

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
# print("===========================================")
# print("==== Be aware: ignoring FutureWarnings ====")
# print("===========================================")

# =================================== FUNC ====================================
def get_time_shift(event_id):
    """get the time shift between centroid time and hypocenter time to
    shift the synthetic seismogram into absolute time using the equation
    t_abs = t_pde + time shift + t_syn
    """
    MT = getdata.get_GCMT_solution(event_id)
    CMTSOLUTIONPATH = (pathnames()['kupedata'] +
                                'CMTSOLUTIONS/{}CMTSOLUTION'.format(event_id))
    CMTSOLUTION = read_events(CMTSOLUTIONPATH)

    CMTSOLUTION_time = CMTSOLUTION[0].origins[0].time
    CENTROID_time = [i.time for i in MT.origins
                                            if i.origin_type == "centroid"][0]

    time_shift = abs(CMTSOLUTION_time - CENTROID_time)

def initial_data_gather(event_id,code,tmin,tmax):
    """gather event information, observation and synthetic traces,
    preprocess all traces accordingly and return one stream object with 6 traces
    """
    # station information
    network,station,location,channel = code.split('.')
    component = channel[-1]

    # event information
    MT = get_data.get_GCMT_solution(event_id)
    moment_tensor = MT.focal_mechanisms[0].moment_tensor
    half_duration = (moment_tensor.source_time_function['duration'])/2

    # grab synthetic data locally
    syntheticdata_path = join(pathnames()['syns'],event_id,'')
    syntheticdata = Stream()
    for c in ["N","E","Z"]:
        syntheticdata_filename = "{n}.{s}.BX{co}.semv.mseed".format(n=network,
                                                                s=station,
                                                                co=c)
        syntheticdata += read(join(syntheticdata_path,syntheticdata_filename))

    # grab observation data
    observationdata,inv,cat = get_data.event_stream(station=station,
                                            channel=channel,
                                            event_id=event_id,
                                            startpad=0,
                                            endpad=350)

    # preprocessing, instrument response, STF convolution (synthetics)
    observationdata_proc = procmod.preprocess(observationdata,
                                                inv=inv,
                                                output='DISP')

    syntheticdata_preproc = procmod.preprocess(syntheticdata)
    syntheticdata_proc = synmod.stf_convolve(syntheticdata_preproc,
                                                            half_duration)

    # combine, common sampling rate, filter, trim common time
    st_IDG = observationdata_proc + syntheticdata_proc
    st_IDG.filter('bandpass',freqmin=1/tmax,
                             freqmax=1/tmin,
                             corners=2,
                             zerophase=True)
    procmod.trimstreams(st_IDG)

    return st_IDG

def plot_obsynth(st,twinax=False):
    """plot 6 streams in 3 subplot figure
    """
    f,(ax1,ax2,ax3) = plt.subplots(3,sharex=True,sharey=False,
                                                    figsize=(9,5),dpi=200)
    # create time axis
    stats = st[0].stats
    t = np.linspace(0,stats.endtime-stats.starttime,stats.npts)

    # axes lists for plotting
    axes = [ax1,ax2,ax3]
    twin_axes = [ax1,ax2,ax3]
    if twinax:
        ax1a,ax2a,ax3a = ax1.twinx(),ax2.twinx(),ax3.twinx()
        twin_axes = [ax1a,ax2a,ax3a]
        ax2a.set_ylabel('BX* velocity (m/s)')

    obs_list = ["HHN","HHE","HHZ"]
    syn_list = ["BXN","BXE","BXZ"]

    # plot
    for ax,tax,obs,syn in zip(axes,twin_axes,obs_list,syn_list):
        obs_select = st.select(channel=obs)[0]
        syn_select = st.select(channel=syn)[0]
        A = ax.plot(t,obs_select.data,color='k',label=obs_select.get_id())
        B = tax.plot(t,syn_select.data,color='r',label=syn_select.get_id())
        plotmod.pretty_grids(ax)
        lines = A+B
        labels = [l.get_label() for l in lines]
        ax.legend(lines,labels,prop={"size":5})
        plotmod.pretty_grids(tax)
        plotmod.align_yaxis(ax,0,tax,0)

    # final plot
    ax1.set_title(event_id)
    ax2.set_ylabel('HH* velocity (m/s)')
    ax3.set_xlabel('time (sec)')
    plt.show()

    return f

def perl_vs_personal(st):
    """plot 6 streams in 3 subplot figure
    """
    f,(ax1,ax2,ax3) = plt.subplots(3,sharex=True,sharey=False,
                                                    figsize=(9,5),dpi=200)
    # create time axis
    stats = st[0].stats
    t = np.linspace(0,stats.endtime-stats.starttime,stats.npts)

    # axes lists for plotting
    axes = [ax1,ax2,ax3]
    obs_list = ["HHN","HHE","HHZ"]
    syn_list = ["BXN","BXE","BXZ"]

    # plot    # check
    for ax,syn in zip(axes,syn_list):
        obs_streams = st.select(location="")
        obs_select = obs_streams.select(channel=syn)[0]
        syn_streams = st.select(location="SAC")
        syn_select = syn_streams.select(channel=syn)[0]
        A = ax.plot(t,obs_select.data,color='k',label=obs_select.get_id())
        B = ax.plot(t,syn_select.data,color='r',label=syn_select.get_id())
        ax.legend(prop={"size":5})
        plotmod.pretty_grids(ax)

    # final plot
    ax1.set_title(event_id)
    ax2.set_ylabel('velocity (m/s)')
    ax3.set_xlabel('time (sec)')
    plt.show()
    stats = st[0].stats

        # import glob
        # st_sac_file = glob.glob('/seis/prj/fwi/bchow/spectral/kupe/test_perl/*.bp')
        # sac_stream = Stream()
        # for sac in st_sac_file:
        #     sac_file = read(sac)
        #     sac_file[0].stats.location = 'SAC'
        #     sac_stream += sac_file
        # st_new = st + sac_stream
        # st_new.resample(50)
        # fig = perl_vs_personal(st_new)

    return f


# =================================== MAIN ====================================
if __name__ == "__main__":
    # station_code = sys.argv[1].upper()
    station_name_path = (pathnames()['data'] +
                                'STATIONXML/north_island_BB_station_names.npy')
    station_names = np.load(station_name_path)
    event_id = "2014p240655"

    for station_code in station_names:
        try:
            code = "NZ.{}..HH?".format(station_code)
            st = initial_data_gather(event_id,code,tmin=5,tmax=30)
            fig = plot_obsynth(st,twinax=False)
        except Exception as e:
            print(e)
            continue

    # import ipdb;ipdb.set_trace()
