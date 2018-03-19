import os
import sys
sys.path.append('../modules/')
import numpy as np
from os.path import join
from obspy import UTCDateTime, read, Stream

# module functions
import getdata
import synmod
import procmod
import plotmod
from getdata import pathnames

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
def initial_data_gather(code,event_id,tmin,tmax):
    """gather event information, observation and synthetic traces,
    preprocess all traces accordingly and return one stream object with 6 traces
    """
    # station information
    net,sta,loc,cha = code.split('.')

    # event information
    time_shift, half_duration = synmod.tshift_halfdur(event_id)

    # grab synthetic data locally
    syntheticdata_path = join(pathnames()['syns'],event_id,'')
    syntheticdata = Stream()
    for c in ["N","E","Z"]:
        syntheticdata_filename = "{n}.{s}.BX{co}.semv.mseed".format(n=net,
                                                                s=sta,
                                                                co=c)
        syntheticdata += read(join(syntheticdata_path,syntheticdata_filename))

    # grab observation data
    observationdata,inv,cat = getdata.event_stream(code=code,
                                                    event_id=event_id,
                                                    startpad=0,
                                                    endpad=350)

    # preprocessing, instrument response, STF convolution (synthetics)
    observationdata_proc = procmod.preprocess(observationdata,
                                                inv=inv,
                                                output='DISP')
    syntheticdata_preproc = synmod.stf_convolve(st=syntheticdata,
                                                 half_duration=half_duration,
                                                 time_shift=time_shift)
    syntheticdata_proc = procmod.preprocess(syntheticdata_preproc)


    # combine, common sampling rate, filter, trim common time
    st_IDG = observationdata_proc + syntheticdata_proc
    procmod.trimstreams(st_IDG)

    st_IDG.filter('bandpass',freqmin=1/tmax,
                             freqmax=1/tmin,
                             corners=2,
                             zerophase=True)

    return st_IDG

def plot_obsynth(st,twinax=False,save=False):
    """plot 6 streams in 3 subplot figure
    """
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

    obs_list = ["HHN","HHE","HHZ"]
    if net == "YH":
        obs_list = ["HH1","HH2","HHZ"]
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
    ax1.set_xlim([t.min(),t.max()])
    ax1.set_title("{e} {s}".format(e=event_id,s=stats.station))
    ax2.set_ylabel("velocity (m/s)")
    ax3.set_xlabel("time (sec)")

    if save:
        figtitle = "{e}_{s}.png".format(e=event_id,s=stats.station)
        figfolder = join(pathnames()['kupeplots'],"obsynth_plots",figtitle)
        plt.savefig(figfolder,dpi=200)

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
    station_names = ["LOBS1","LOBS2","LOBS3","LOBS4"]
    event_id = "2014p715167"

    for station_code in station_names:
        # try:
        # code = "NZ.{}.*.HH?".format(station_code)
        code = "YH.{}.*.HH?".format(station_code)
        st = initial_data_gather(code,event_id,tmin=5,tmax=30)
        fig = plot_obsynth(st,twinax=True,save=False)
        # except Exception as e:
        #     print(e)
        #     continue
