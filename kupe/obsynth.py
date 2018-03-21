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

def plot_event_station(inv,cat):
    """plot event and station along with line connection and station-receiver
    distance for quick visualization
    """

    # plot station on a map
    fig = inv.plot(marker='v',
            projection="local",
            resolution = "i",
            size=50,
            show=False,
            continent_fill_color="white",
            water_fill_color="white",
            color_per_network=True,
            label=False)

    # convert station and event information to map coords
    station_lat = inv[0][0].latitude
    station_lon = inv[0][0].longitude
    station_x,station_y = fig.bmap(station_lon,station_lat)
    event_lat = cat[0].origins[0].latitude
    event_lon = cat[0].origins[0].longitude
    event_x,event_y = fig.bmap(event_lon,event_lat)

    # annotate station code
    stationcode = inv[0][0].code
    plt.annotate(stationcode,
                xy=(station_x,station_y),
                xytext=(station_x,station_y),
                fontsize=10,
                weight='bold',
                zorder=100)

    # plot event
    scatter = fig.bmap.scatter(event_x,event_y,
                                marker='o',
                                s=150,
                                zorder=50,
                                edgecolor='k')
    magnitude = round(cat[0].magnitudes[0].mag,2)
    plt.annotate(magnitude,
                xy=(event_x,event_y),
                xytext=(event_x,event_y),
                fontsize=10,
                weight='bold',
                zorder=100)

    # connect station and event by arrow
    dist_azi = gps2dist_azimuth(event_lat,event_lon,station_lat,station_lon)
    epi_dist = round(dist_azi[0],2)
    BAz = round(dist_azi[2],2)
    plt.title("Epicentral Distance: {} | BAz: {}".format(epi_dist,BAz))


    plt.show()


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
    # ======================== GCMT ADDITION ==========================
    syntheticdata_path_PART2 = join(pathnames()['syns'],event_id,
                                                    "with_GCMT_solution",'')
    # ======================== GCMT ADDITION ==========================

    syntheticdata = Stream()
    for c in ["N","E","Z"]:
        syntheticdata_filename = "{n}.{s}.BX{co}.semv.mseed".format(n=net,
                                                                s=sta,
                                                                co=c)
        syntheticdata += read(join(syntheticdata_path,syntheticdata_filename))
        # ======================== GCMT ADDITION ==========================
        syntheticdata_PART2 = read(join(syntheticdata_path_PART2,
                                                        syntheticdata_filename))
        syntheticdata_PART2[0].stats.location = "withGCMT"
        syntheticdata += syntheticdata_PART2
        # ======================== GCMT ADDITION ==========================

    # grab observation data
    observationdata,inv,cat = getdata.event_stream(code=code,
                                                    event_id=event_id,
                                                    startpad=0,
                                                    endpad=350)

    # plot event and station on a map
    if plotmap:
        plot_event_station(inv,cat)


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
    if twinax:
        ax1a,ax2a,ax3a = ax1.twinx(),ax2.twinx(),ax3.twinx()
        twin_axes = [ax1a,ax2a,ax3a]
        ax2a.set_ylabel("velocity (m/s)")

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
        # syn_select = st.select(channel=syn)[0]
        # ======================== GCMT ADDITION ==========================
        syn_select = st.select(channel=syn)
        syn_select_PART2 = syn_select.select(location="withGCMT")[0]
        syn_select = syn_select.select(location="")[0]
        C = tax.plot(t,syn_select_PART2.data,color='r',label=syn_select_PART2.get_id(),
                                                linestyle='--')
        # ======================== GCMT ADDITION ==========================

        A = ax.plot(t,obs_select.data,color='k',label=obs_select.get_id())
        B = tax.plot(t,syn_select.data,color='r',label=syn_select.get_id())
        plotmod.pretty_grids(ax)
        lines = A+B+C
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
if __name__ == "__main__":
    # station_code = sys.argv[1].upper()
    station_name_path = (pathnames()['data'] +
                                'STATIONXML/NAMESOF_NZ_NI_BB_seismometers.npy')
    station_names = np.load(station_name_path)
    station_names = ['BFZ','BKZ','HAZ','HIZ','KNZ','MRZ','MWZ','OPRZ','PUZ',
                        'PXZ','RTZ','TLZ','TOZ','TSZ','VRZ','WAZ']
    event_id = "2014p240655"

    for station_code in station_names:
        try:
            code = "NZ.{}.*.HH?".format(station_code)
            bounds = [6,30]
            st = initial_data_gather(code,event_id,
                                    bounds=bounds,
                                    plotmap=False)
            fig = plot_obsynth(st,bounds=bounds,
                                    twinax=False,
                                    save=True,
                                    show=False)
        except Exception as e:
            print("="*79)
            traceback.print_exc()
            print("="*79)
            continue
