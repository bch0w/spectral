"""script for building adjoint source using the python packages pyflex
and pyadjoint. includes preprocessing of observations and synthetics so that
they are in the proper configuration for feeding into aformentioned packages.
+initial code for data gather taken from obsynth.py
"""
import os
import sys
import pprint
import pyflex
import pyadjoint
import numpy as np
from os.path import join
from obspy import UTCDateTime, read, Stream

# module functions
sys.path.append('../modules/')
import getdata
import synmod
import procmod
import plotmod
from getdata import pathnames

# plotting settings
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=mpl.cbook.mplDeprecation)

# ============================ HELPER FUNCTIONS ================================
def find_BAz(inv,event):
    """get backazimuth based on obspy inventory and event information
    """
    from obspy.geodetics import gps2dist_azimuth
    station_lat = inv[0][0].latitude
    station_lon = inv[0][0].longitude
    event_lat = event.origins[0].latitude
    event_lon = event.origins[0].longitude
    dist,Az,BAz = gps2dist_azimuth(event_lat,event_lon,station_lat,station_lon)

    return BAz

def breakout_stream(st):
    """get observed and synthetic parts of a stream object
    assumes that 'st' contains all observation and synthetic traces and that
    they follow the naming scheme: obs="HH?"", syn="BX?"

    :type st: obspy.stream
    :param st: stream to be separated
    :rtype obs_stream: obspy.stream
    :return obs_stream: single trace stream containing observation data
    :rtype syn_stream: obspy.stream
    :return syn_stream: single trace stream contianing synthetic data
    """
    # break out streams into individuals
    obs_stream = st.select(channel="HH?")
    syn_stream = st.select(channel="BX?")

    return obs_stream, syn_stream

def get_station_latlon(sta):
    """geonet response files dont have station information, so grab it internal
    """
    # find station coordinates
    nz_bb_path = pathnames()['data'] + 'STATIONXML/nz_BB_coords.npz'
    nz_bb = np.load(nz_bb_path)
    nz_bb_names = nz_bb['NAME']

    bb_ind = np.where(nz_bb_names==sta)[0][0]
    sta_lat = nz_bb['LAT'][bb_ind]
    sta_lon = nz_bb['LON'][bb_ind]

    return sta_lat, sta_lon

# ============================= MAIN FUNCTIONS =================================
def initial_data_gather(PD):
    """gather event information, observation and synthetic traces.
    preprocess all traces accordingly and return one stream object with 6 traces
    stolen and modified from obsynth.py

    TODO:
        +by default grabs all components, could be modified to only grab
        necessary components, if the speedup is worthwhile
        +change the stream length by the theoretical arrival times or
        source-receiver distance (?)


    :type PD: dictionary
    :param PD: parameter dictionary
    :rtype st_IDG: obspy.stream
    :return st_IDG: fully preprocessed observed and synthetic data, ready to go
    :rtype inv: obspy.inventory
    :return inv: response and station information
    :rtype event: obspy.event
    :return event: event information
    """
    # station information
    net,sta,loc,cha = PD["code"].split('.')

    # grab synthetic data locally, decide
    syntheticdata_path = join(pathnames()['syns'],PD["event_id"],'')

    syntheticdata = Stream()
    for c in ["N","E","Z"]:
        syntheticdata_filename = "{n}.{s}.BX{co}.semv.mseed".format(n=net,
                                                                    s=sta,
                                                                    co=c)
        syntheticdata += read(join(syntheticdata_path,syntheticdata_filename))

    # grab observation data
    observationdata,inv,cat = getdata.event_stream(code=PD["code"],
                                                    event_id=PD["event_id"],
                                                    startpad=0,
                                                    endpad=180)
    event = cat[0]
    if inv[0][0].latitude == 0.0:
        sta_lat,sta_lon = get_station_latlon(sta)
        inv[0][0].latitude = sta_lat
        inv[0][0].longitude = sta_lon

    if not observationdata:
        print("No observation data")
        return None,None,None

    # rotate to theoretical backazimuth if necessary
    if PD["component"] == ("R" or "T"):
        BAz = find_BAz(inv,event)
        observationdata_proc.rotate(method='NE->RT',back_azimuth=BAz)
        syntheticdata_proc.rotate(method='NE->RT',back_azimuth=BAz)

    # only grab the necessary streams to save on processing
    observationdata = observationdata.select(component=PD["component"])
    syntheticdata = syntheticdata.select(component=PD["component"])

    # observation preprocessing + instrument response
    observationdata = procmod.preprocess(observationdata,
                                                inv=inv,
                                                output=PD["output"])

    # synthetic moment-tensor information
    time_shift, half_duration = synmod.tshift_halfdur(PD["event_id"])

    # if GCMT solution doesn't exist, timeshift isn't possible
    if time_shift:
        syntheticdata = synmod.stf_convolve(st=syntheticdata,
                                                 half_duration=half_duration,
                                                 time_shift=time_shift)

    syntheticdata = procmod.preprocess(syntheticdata,inv=None,
                                                        output=PD["output"])


    # combine and trim to common time
    st_IDG = observationdata + syntheticdata
    st_IDG = procmod.trimstreams(st_IDG)

    # filter
    tmin,tmax = PD["bounds"]
    st_IDG.filter('bandpass',freqmin=1/tmax,
                             freqmax=1/tmin,
                             corners=2,
                             zerophase=True)

    return st_IDG,inv,event

def run_pyflex(PD,st,inv,event,plot=False,config="UAF"):
    """use pyflex to grab windows, current config set to defaults found on docs

    :type PD: dictionary
    :param PD: parameter sets
    :rtype windows: pyflex.window object
    :return windows: windows containing selected timespans where waveforms
    have acceptable match, also contains information about the window (see docs)
    """
    obs,syn = breakout_stream(st)

    # DEFAULT CONFIG
    if config == "default":
        config = pyflex.Config(min_period=PD["bounds"][0],
                               max_period=PD["bounds"][1],
                               stalta_waterlevel=0.08,
                               tshift_acceptance_level=15.0,
                               dlna_acceptance_level=1.0,
                               cc_acceptance_level=0.8,
                               c_0=0.7,
                               c_1=4.0,
                               c_2=0.0,
                               c_3a=1.0,
                               c_3b=2.0,
                               c_4a=3.0,
                               c_4b=10.0
                               )
        # UAF CONFIG
    elif config == "UAF":
        config = pyflex.Config(min_period=PD["bounds"][0],
                               max_period=PD["bounds"][1],
                               stalta_waterlevel=0.18,
                               tshift_acceptance_level=4.0,
                               dlna_acceptance_level=1.5,
                               cc_acceptance_level=0.71,
                               c_0=0.7,
                               c_1=2.0,
                               c_2=0.0,
                               c_3a=3.0,
                               c_3b=2.0,
                               c_4a=2.5,
                               c_4b=12.0)

    pf_event = pyflex.Event(latitude=event.origins[0].latitude,
                            longitude=event.origins[0].longitude,
                            depth_in_m=event.origins[0].depth,
                            origin_time=event.origins[0].time)

    pf_station = pyflex.Station(latitude=inv[0][0].latitude,
                                longitude=inv[0][0].longitude)

    windows = pyflex.select_windows(observed=obs,
                                    synthetic=syn,
                                    config=config,
                                    event=pf_event,
                                    station=pf_station,
                                    plot=plot)
    if not windows:
        print("Empty windows")
        return None

    return windows

# def pyflex_window_viewer(obs,syn,windows):


def run_pyadjoint(PD,st,windows,output_path=None,plot=False):
    """function to call pyadjoint with preset configurations
    ++not in the docs:
    in pyadjoint.calculate_adjoint_source: window needs to be a list of lists,
    with each list containing the [left_window,right_window] where each window
    argument is given in seconds

    :type PD: dictionary
    :param PD: parameter sets
    :type st: obspy.stream
    :param st: stream to be separated
    :type windows: pyflex.window object
    :param windows: windows containing selected timespans where waveforms
    have acceptable match, also contains information about the window (see docs)
    :type output_path: str
    :param output_path: where to save the adjoint source, if none, no saving
    :rtpye adj_src: pyadjoint
    """
    obs,syn = breakout_stream(st)
    net,sta,loc,cha = PD["code"].split('.')
    cha=cha[:2]+PD["component"]

    # collect all windows into a single list object
    adjoint_windows = []
    delta = st[0].stats.delta
    for win in windows:
        adj_win = [win.left*delta,win.right*delta]
        adjoint_windows.append(adj_win)

    config = pyadjoint.ConfigWaveForm(min_period=PD["bounds"][0],
                                      max_period=PD["bounds"][1],
                                      taper_type="hann",
                                      taper_percentage=0.15)
    config_adj = pyadjoint.ConfigCrossCorrelation(min_period=PD["bounds"][0],
                                                  max_period=PD["bounds"][1],
                                                  taper_type='hann',
                                                  taper_percentage=0.3,
                                                  measure_type='dt',
                                                  dt_sigma_min=1.0,
                                                  dlna_sigma_min=0.5)


    adj_src = pyadjoint.calculate_adjoint_source(
                                 adj_src_type="cc_traveltime_misfit",
                                 observed=obs,
                                 synthetic=syn,
                                 config=config_adj,
                                 window=adjoint_windows,
                                 plot=plot
                                 )

    if output_path:
        output_file = os.path.join(output_path,
                                   "{NET}.{STA}.{CHA}.{EVENT_ID}.adj".format(
                                                     NET=net,
                                                     STA=sta,
                                                     CHA=cha,
                                                     EVENT_ID=PD["event_id"]))
        adj_src.write(output_file,
                      format="SPECFEM",
                      time_offset=0)

    return adj_src

def bobTheBuilder():
    """main processing script

    ++intra-function parameters and choices:
    :type EVENT_ID: str
    :param EVENT_ID: GEONET event ID
    :type STATION_NAME: str
    :param STATION_NAME: identifier for data fetching, in the form NN.SSS(S)
                         examples: NZ.HIZ, XX.RD01
    :type MINIMUM/MAXIMUM_FILTER_PERIOD: int
    :param MINIMUM/MAXIMUM_FILTER_PERIOD: bandpass filtering bounds in seconds
    :type COMPONENT: str
    :param COMPONENT: component of choice to be used, available N,E,Z,T,R
    :type UNIT_OUTPUT: str
    :param UNIT_OUTPUT: available DISP, VEL, ACC
    :type ADJ_SRC_OUTPUT_PATH: str
    :param ADJ_SRC_OUTPUT_PATH: "path/to/save/adjoint_source/"
    :type PLOT: bool
    :param PLOT: plot outputs of pyflex and pyadjoint, global switch
    """
    # =============== PARAMETER SET ===============
    EVENT_IDS = ["2018p130600"]
    STATION_NAMES = ['NZ.BFZ','NZ.BKZ','NZ.HAZ','NZ.HIZ','NZ.KNZ','NZ.MRZ',
                     'NZ.MWZ','NZ.OPRZ','NZ.PUZ','NZ.PXZ','NZ.RTZ','NZ.TLZ',
                     'NZ.TOZ','NZ.TSZ','NZ.VRZ','NZ.WAZ']
    # STATION_NAMES = ['XX.RD01','XX.RD02','XX.RD03','XX.RD04','XX.RD05',
    #                    'XX.RD06','XX.RD07','XX.RD08','XX.RD09','XX.RD10',
    #                    'XX.RD11','XX.RD12','XX.RD13','XX.RD14','XX.RD15',
    #                    'XX.RD16','XX.RD17','XX.RD18','XX.RD19','XX.RD20',
    #                    'XX.RD21','XX.RD22']
    MINIMUM_FILTER_PERIOD = 6
    MAXIMUM_FILTER_PERIOD = 30
    COMPONENT = "Z"
    UNIT_OUTPUT = "VEL"
    # ADJOINT_TYPE = "cc_traveltime_misfit"
    # CONFIG = "default"
    ADJ_SRC_OUTPUT_PATH = pathnames()["kupedata"] + "ADJOINTSOURCES"
    ADJ_SRC_OUTPUT_PATH = None
    PLOT = True
    # =============== PARAMETER SET ===============

    for EVENT_ID in EVENT_IDS:
        for STATION_NAME in STATION_NAMES:
            print(STATION_NAME)
            PAR_DICT = {"station_name":STATION_NAME,
                        "code":"{}.*.HH?".format(STATION_NAME),
                        "event_id":EVENT_ID,
                        "bounds":(MINIMUM_FILTER_PERIOD,
                                  MAXIMUM_FILTER_PERIOD),
                        "component":COMPONENT,
                        "output":UNIT_OUTPUT
                        }

            # PROCESSING
            st,inv,event = initial_data_gather(PAR_DICT)
            if not st:
                continue
            windows = run_pyflex(PAR_DICT,st,inv,event,plot=PLOT)
            import ipdb;ipdb.set_trace()
            if not windows:
                continue
            adj_src = run_pyadjoint(PAR_DICT,st,windows,
                                    output_path=ADJ_SRC_OUTPUT_PATH,
                                    plot=PLOT)



# =================================== MAIN ====================================
if __name__ == "__main__":
    bobTheBuilder()
