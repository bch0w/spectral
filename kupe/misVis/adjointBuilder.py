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

import windowMaker
import mapMaker

# module functions
sys.path.append('../../modules/')
import getdata
import synmod
import procmod
import plotmod
from getdata import pathnames

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

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
    nz_bb_path = pathnames()['data'] + 'STATIONXML/NZ_BB_coords.npz'
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
                                                    startpad=20,
                                                    endpad=180)
    event = cat[0]
    if inv[0][0].latitude == 0.0:
        sta_lat,sta_lon = get_station_latlon(sta)
        inv[0][0].latitude = sta_lat
        inv[0][0].longitude = sta_lon

    if not observationdata:
        print("No observation data")
        return None,None,None

    # observation preprocessing + instrument response
    observationdata = procmod.preprocess(observationdata,
                                                inv=inv,
                                                output=PD["output"])
    
    # rotate to theoretical backazimuth if necessary
    if PD["rotate"] == True:
        BAz = find_BAz(inv,event)
        observationdata.rotate(method='NE->RT',back_azimuth=BAz)
        syntheticdata.rotate(method='NE->RT',back_azimuth=BAz)

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

def choose_config(config):
    """helper function to avoid typing out the full pyflex config, stores values
    in a list with the following order
    0:stalta_Waterlevel, 1:tshift_acceptance_level, 2:dlna_acceptance_level,
    3:cc_acceptance_level, 4:c_0, 5:c_1, 6:c_2, 7:c_3a, 8:c_3b, 10:c_4a, 11:c_4b
    """
    cfgdict = {"default":[.08,15.,1.,.8,.7,4.,0.,1.,2.,3.,10.],
                  "UAF":[.18,4.,1.5,.71,.7,2.,0.,3.,2.,2.5,12.]}

    return cfgdict[config]

def run_pyflex(PD,st,inv,event,plot=False):
    """use pyflex to grab windows, current config set to defaults found on docs

    :type PD: dictionary
    :param PD: parameter sets
    :rtype windows: pyflex.window object
    :return windows: windows containing selected timespans where waveforms
    have acceptable match, also contains information about the window (see docs)
    """
    CD = choose_config(PD["pyflex_config"])
    config = pyflex.Config(min_period=PD["bounds"][0],
                           max_period=PD["bounds"][1],
                           stalta_waterlevel=CD[0],
                           tshift_acceptance_level=CD[1],
                           dlna_acceptance_level=CD[2],
                           cc_acceptance_level=CD[3],
                           c_0=CD[4],c_1=CD[5],c_2=CD[6],c_3a=CD[7],
                           c_3b=CD[8],c_4a=CD[9],c_4b=CD[10])

    pf_event = pyflex.Event(latitude=event.origins[0].latitude,
                            longitude=event.origins[0].longitude,
                            depth_in_m=event.origins[0].depth,
                            origin_time=event.origins[0].time)

    pf_station = pyflex.Station(latitude=inv[0][0].latitude,
                                longitude=inv[0][0].longitude)

    # iterate windows by component and place into dictionary output
    windows = {}
    for comp in PD["comp_list"]:
        obs,syn = breakout_stream(st.select(component=comp))
        window = pyflex.select_windows(observed=obs,
                                        synthetic=syn,
                                        config=config,
                                        event=pf_event,
                                        station=pf_station,
                                        plot=plot)
        # stalta = pyflex.stalta.sta_lta(data=syn[0].data,
        #                                dt=syn[0].stats.delta,
        #                                min_period=PD["bounds"][0])
        windows[comp] = window
    
    # staltas['waterlevel'] = CD[0]      
                              
    if not windows:
        print("Empty windows")
        return None

    return windows


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

def build_figure(st,inv,event,windows,PD):
    """take outputs of mapMaker and windowMaker and put them into one figure
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    
    # attempting to plot everything on one figure, mpl doesn't like
    # F = plt.figure(figsize=(11.69,8.27),dpi=100)
    # width_ratios = [1,1]
    # GS = mpl.gridspec.GridSpec(1,2,width_ratios=width_ratios)
    # 
    # F1 = F.add_subplot(GS[0])
    # axes = windowMaker.window_maker(st,windows,PD=PD)
    # 
    # F2 = F.add_subplot(GS[1])
    # axes = windowMaker.window_maker(st,windows,PD=PD)
    
    f2 = plt.figure(figsize=(11.69,8.27),dpi=100)
    axes = windowMaker.window_maker(st,windows,PD=PD)
    
    f2 = plt.figure(figsize=(10,9.4),dpi=100)
    map = mapMaker.generate_map(event,inv,faults=PD['plot_faults'])
    
    plt.show()
    
def bob_the_builder():
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
    # ============================== PARAMETER SET =============================
    
    EVENT_IDS = ["2014p240655"]
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
    ROTATE = True
    UNIT_OUTPUT = "VEL"
    PYFLEX_CONFIG = "UAF"
    ADJ_SRC_OUTPUT_PATH = pathnames()["kupedata"] + "ADJOINTSOURCES"
    PLOT = False
    PLOT_FAULTS = True
    SAVE_PLOT = True
    # ADJOINT_TYPE = "cc_traveltime_misfit"

    # ============================== PARAMETER SET =============================
    
    # PARAMETER DEFUALT SET
    COMPONENT_LIST = ["N","E","Z"]
    if ROTATE:
        COMPONENT_LIST = ["R","T","Z"]
    
    # MAIN ITERATE OVER EVENTS
    for EVENT_ID in EVENT_IDS:
        for STATION_NAME in STATION_NAMES:
            print(STATION_NAME)
            PAR_DICT = {"station_name":STATION_NAME,
                        "code":"{}.*.HH?".format(STATION_NAME),
                        "event_id":EVENT_ID,
                        "bounds":(MINIMUM_FILTER_PERIOD,
                                  MAXIMUM_FILTER_PERIOD),
                        "rotate":ROTATE,
                        "output":UNIT_OUTPUT,
                        "pyflex_config":PYFLEX_CONFIG,
                        "comp_list":COMPONENT_LIST,
                        "save_plot":SAVE_PLOT,
                        "plot_faults":PLOT_FAULTS
                        }
            
            # MAIN PROCESSING
            st,inv,event = initial_data_gather(PAR_DICT)
            if not st: continue
            windows = run_pyflex(PAR_DICT,st,inv,event,plot=PLOT)
            if not windows: continue
            build_figure(st,inv,event,windows,PAR_DICT)
            # adj_src = run_pyadjoint(PAR_DICT,st,windows,
            #                         output_path=ADJ_SRC_OUTPUT_PATH,
            #                         plot=PLOT)
                                                                     

def _test_build_figure():
    """test figure building with example data
    """
    from obspy import read_events, read_inventory
    boundsdict = {"station_name":"TEST","bounds":(6,30)}
    streampath = pathnames()['data'] + 'WINDOWTESTING/testmseed.pickle'
    windowpath = pathnames()['data'] + 'WINDOWTESTING/testwindows.npz'
    st = read(streampath)
    windows = np.load(windowpath)
    eventpath = pathnames()['data'] + "WINDOWTESTING/testevent.xml"
    invpath = pathnames()['data'] + "WINDOWTESTING/testinv.xml"
    cat = read_events(eventpath)
    event = cat[0]
    inv = read_inventory(invpath)
    
    build_figure(st,inv,event,windows,boundsdict)

# =================================== MAIN ====================================
if __name__ == "__main__":
    # _test_build_figure()
    bob_the_builder()

