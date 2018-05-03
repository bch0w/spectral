import os
import sys
import traceback
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
warnings.filterwarnings("ignore", category=mpl.cbook.mplDeprecation)

# print("===========================================")
# print("==== Be aware: ignoring FutureWarnings ====")
# print("===========================================")

# =================================== FUNC ====================================
def find_BAz(inv,cat):
    """get backazimuth
    """
    station_lat = inv[0][0].latitude
    station_lon = inv[0][0].longitude
    event_lat = cat[0].origins[0].latitude
    event_lon = cat[0].origins[0].longitude
    BAz = gps2dist_azimuth(event_lat,event_lon,station_lat,station_lon)[2]

    return BAz


def initial_data_gather(code,event_id,bounds,plotmap=False):
    """gather event information, observation and synthetic traces,
    preprocess all traces accordingly and return one stream object with 6 traces
    """
    # station information
    net,sta,loc,cha = code.split('.')

    # event information
    time_shift, half_duration = synmod.tshift_halfdur(event_id)

    # filter bounds
    tmin,tmax = bounds

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
    # plot event and station on a map
    if plotmap:
        plotmod.plot_event_station(inv,cat)


    # preprocessing, instrument response, STF convolution (synthetics)
    observationdata_proc = procmod.preprocess(observationdata,
                                                inv=inv,
                                                output='DISP')
    syntheticdata_preproc = synmod.stf_convolve(st=syntheticdata,
                                                 half_duration=half_duration,
                                                 time_shift=time_shift)
    syntheticdata_proc = procmod.preprocess(syntheticdata_preproc)

    # rotate to theoretical backazimuth
    # BAz = find_BAz(inv,cat)
    # observationdata_proc.rotate(method='NE->RT',back_azimuth=BAz)
    # syntheticdata_proc.rotate(method='NE->RT',back_azimuth=BAz)

    # combine, common sampling rate, filter, trim common time
    st_IDG = observationdata_proc + syntheticdata_proc
    procmod.trimstreams(st_IDG)

    st_IDG.filter('bandpass',freqmin=1/tmax,
                             freqmax=1/tmin,
                             corners=2,
                             zerophase=True)

    return st_IDG

def plot_obsynth(st,bounds,twinax=False,save=False,show=True):
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
    ax2.set_ylabel("velocity (m/s)\n{}".format(c2))
    if twinax:
        ax1a,ax2a,ax3a = ax1.twinx(),ax2.twinx(),ax3.twinx()
        twin_axes = [ax1a,ax2a,ax3a]
        ax2a.set_ylabel("synthetic velocity (m/s)")
        ax2.set_ylabel("observed velocity (m/s)\n{}".format(c2))



    obs_list = ["HH{}".format(c1),"HH{}".format(c2),"HH{}".format(c3)]
    if net == "YH":
        obs_list = ["HHZ","HH1","HH2"]
    syn_list = ["BX{}".format(c1),"BX{}".format(c2),"BX{}".format(c3)]
    # obs_list = ["HHZ","HHR","HHT"]
    # if net == "YH":
    #     obs_list = ["HH1","HH2","HHZ"]
    # syn_list = ["BXZ","BXR","BXT"]

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
    ax1.set_title("{e} {s} [{t0}-{t1}s]".format(e=event_id,
                                                    s=stats.station,
                                                    t0=bounds[0],
                                                    t1=bounds[1])
                                                    )

    ax1.set_ylabel(c1)
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
if __name__ == "__main__":
    # station_code = sys.argv[1].upper()
    station_name_path = (pathnames()['data'] +
                                'STATIONXML/NAMESOF_NZ_NI_BB_seismometers.npy')
    station_names = np.load(station_name_path)
    station_names = ['BFZ','BKZ','HAZ','HIZ','KNZ','MRZ','MWZ','OPRZ','PUZ',
                        'PXZ','RTZ','TLZ','TOZ','TSZ','VRZ','WAZ']

    event_id = "2016p669820"

    for station_code in station_names:
        try:
            code = "NZ.{}.*.HH?".format(station_code)
            bounds = [6,30]
            st = initial_data_gather(code,event_id,
                                    bounds=bounds,
                                    plotmap=False)
            fig = plot_obsynth(st,bounds=bounds,
                                    twinax=False,
                                    save=False,
                                    show=True)
        except Exception as e:
            print("="*79)
            traceback.print_exc()
            print("="*79)
            continue
