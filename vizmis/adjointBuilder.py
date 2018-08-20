"""script for building adjoint source using the python packages pyflex
and pyadjoint. includes preprocessing of observations and synthetics so that
they are in the proper configuration for feeding into aformentioned packages.
+initial code for data gather taken from obsynth.py
"""
import os
import sys
import pprint
import pyflex
import pyasdf
import pyadjoint
import numpy as np
from os.path import join
from obspy import UTCDateTime, read, Stream
from obspy.signal.filter import envelope

# internal scripts
import windowMaker
import mapMaker
sys.path.append('./tests')
from tests import func_test

# module functions
sys.path.append('../modules/')
import getdata
import synmod
import procmod
import plotmod
from getdata import pathnames

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# ============================ HELPER FUNCTIONS ================================
def _create_log_fid(path):
    """create log identifiers sequentially based on whats available
    """
    import glob
    logtemplate = "LOG_{:0>3}"
    files = glob.glob(join(path,'*'))
    if files:
        lastlognumber = int(max(files).split('_')[1])
    else:
        lastlognumber = -1
    fidout = logtemplate.format(lastlognumber+1)
    
    return fidout

def _check_path(path):
    """small function to check if a path exists and if not, then to create it
    """
    if not os.path.exists(path):
        os.makedirs(path)
    

def _find_BAz(inv,event):
    """get backazimuth based on obspy inventory and event information
    """
    from obspy.geodetics import gps2dist_azimuth
    station_lat = inv[0][0].latitude
    station_lon = inv[0][0].longitude
    event_lat = event.origins[0].latitude
    event_lon = event.origins[0].longitude
    dist,Az,BAz = gps2dist_azimuth(event_lat,event_lon,station_lat,station_lon)

    return BAz

def _get_station_latlon(sta):
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
    
def create_window_dictionary(window):
    """HDF5 doesnt play nice with nonstandard objects in dictionaries, e.g.
    nested dictionaries, UTCDateTime objects. So remake the pyflex window
    json dictionary into something that will sit well in a pyasdf object
    """
    winnDixie = window._get_json_content()

    # change UTCDateTime objects into strings
    winnDixie['absolute_endtime'] = str(winnDixie['absolute_endtime'])
    winnDixie['absolute_starttime'] = str(winnDixie['absolute_starttime'])
    winnDixie['time_of_first_sample'] = str(winnDixie['time_of_first_sample'])

    phase_arrivals = winnDixie['phase_arrivals']
    for phase in phase_arrivals:
        winnDixie['phase_arrival_{}'.format(phase['name'])] = phase['time']

    winnDixie.pop('phase_arrivals')

    return winnDixie

# ============================= MAIN FUNCTIONS =================================
def build_figure(st,inv,event,windows,staltas,adj_src,PD):
    """take outputs of mapMaker and windowMaker and put them into one figure
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    if PD['save_plot'][0]:
        outpath = join(PD['save_plot'][1],PD['event_id'])
        _check_path(outpath)
        
    if PD["plot"][0]:
        if PD["verbose"]: print("Generating waveform plot")
        f1 = plt.figure(figsize=(11.69,8.27),dpi=100)
        axes = windowMaker.window_maker(st,windows,staltas,adj_src,PD=PD)
        if PD['save_plot'][0]:
            if PD["verbose"]:print("Saving figure")
            f1.savefig(join(outpath,'{sta}_wav.png'.format(sta=PD["station"])))

    if PD["plot"][1]:
        if PD["verbose"]: print("Generating source receiver map")
        f2 = plt.figure(figsize=(10,9.4),dpi=100)
        map = mapMaker.generate_map(event,inv,faults=PD["plot"][2])
        if PD['save_plot'][0]:
            if PD["verbose"]:print("Saving figure")
            f2.savefig(join(outpath,'{sta}_map.png'.format(sta=PD["station"])))

    # !!! code snippet to save all open figures into a single PDF, acrobat
    # !!! was having trouble opening up the .pdf maps, but preview can
    # !!! .pdf's much larger than .png so only for profesional outputs
    # import matplotlib.backends.backend_pdf as backend
    # outfid = join(PD['save_plot'][1],
    #               '{id}_{sta}_wavmap.pdf'.format(id=PD["event_id"],
    #                                              sta=PD["station"])
    #                                              )
    # pdf = backend.PdfPages(outfid)
    # for fig in range(1,plt.gcf().number+1):
    #     pdf.savefig(fig)
    # pdf.close()
    
    if PD["show"]:
        plt.show()

    plt.close("all")
    
def initial_data_gather(PD):
    """gather event information, observation and synthetic traces.
    preprocess all traces accordingly and return one stream object with 6 traces
    stolen and modified from obsynth.py

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

    # grab synthetic data locally, by default I was saving mseeds in velocity
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
        sta_lat,sta_lon = _get_station_latlon(sta)
        inv[0][0].latitude = sta_lat
        inv[0][0].longitude = sta_lon

    if not observationdata:
        raise Exception("No observation data")

    # observation preprocessing + instrument response
    if PD["verbose"]:print("Preprocessing observation data")
    observationdata = procmod.preprocess(observationdata,
                                         inv=inv,
                                         output=PD["units"])

    # rotate to theoretical backazimuth if necessary
    if PD["rotate"] == True:
        BAz = _find_BAz(inv,event)
        observationdata.rotate(method='NE->RT',back_azimuth=BAz)
        syntheticdata.rotate(method='NE->RT',back_azimuth=BAz)

    # synthetic moment-tensor information
    time_shift, half_duration = synmod.tshift_halfdur(PD["event_id"])

    # if GCMT solution doesn't exist, timeshift isn't possible
    if time_shift:
        syntheticdata = synmod.stf_convolve(st=syntheticdata,
                                                 half_duration=half_duration,
                                                 time_shift=time_shift)

    if PD["verbose"]:print("Preprocessing synthetic data")
    syntheticdata = procmod.preprocess(syntheticdata,inv=None,
                                       output=PD["units"])


    # combine and trim to common time
    st_IDG = observationdata + syntheticdata
    st_IDG = procmod.trimstreams(st_IDG)

    # filter
    tmin,tmax = PD["bounds"]
    if PD["verbose"]:print("Filtering at {0} to {1} seconds".format(tmin,tmax))
    st_IDG.filter('bandpass',freqmin=1/tmax,
                             freqmax=1/tmin,
                             corners=2,
                             zerophase=True)

    # save into pyasdf dataset if applicable. 'add' function auto writes to file
    # try except statements incase these objects already exist - hacky but
    # excepts catch the errors when trying to add things that are already there
    # add_waveforms will push 'already_exists' exceptions
    if PD["dataset"]:
        if PD["verbose"]:
            print("Saving station, event and waveforms to PyASDF dataset")
        try:
            PD["dataset"].add_quakeml(event)
        except ValueError:
            print('Event already added - exception passed')
            pass
        try:
            PD["dataset"].add_stationxml(inv)
        except TypeError:
            print('Station already added - exception passed')
            pass

        obsout,synout = breakout_stream(st_IDG)
        # ignore ASDFwarnings for already saved data when saving waveforms
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            PD["dataset"].add_waveforms(waveform=obsout,
                                        tag="observed_processed",
                                        event_id=event)
            PD["dataset"].add_waveforms(waveform=synout,
                                        tag="synthetic_processed_{}".format(
                                                                   PD['model']),
                                        event_id=event)

    return st_IDG,inv,event

def choose_config(choice,PD):
    """helper function to avoid typing out the full pyflex or pyadjoint config:
    
    ++PYFLEX (Maggi et al. 2009)
    i  Standard Tuning Parameters:
    0: water level for STA/LTA (short term average/long term average)
    1: time lag acceptance level
    2: amplitude ratio acceptance level (dlna)
    3: normalized cross correlation acceptance level
    i  Fine Tuning Parameters
    4: c_0 = for rejection of internal minima 
    5: c_1 = for rejection of short windows
    6: c_2 = for rejection of un-prominent windows
    7: c_3a = for rejection of multiple distinct arrivals
    8: c_3b = for rejection of multiple distinct arrivals
    9: c_4a = for curtailing windows w/ emergent starts and/or codas
    10:c_4b = for curtailing windows w/ emergent starts and/or codas
    
    ++PYADJOINT:
    different misfit measures detailed in pyadjoint docs, multitaper approach is
    what was used in Tape et al. (2010). All inputs taken from what is specified
    in source code, many are left as default
    """
    if choice == "pyflex":
        cfgdict = {"default":[.08,15.,1.,.8,.7,4.,0.,1.,2.,3.,10.],
                   "UAF":[.18,4.,1.5,.71,.7,2.,0.,3.,2.,2.5,12.],
                   "NZ":[]}
        cfgout = cfgdict[PD["pyflex_config"]]

    elif choice == "pyadjoint":
        if PD["adj_src_type"] == "waveform":
            cfgout = pyadjoint.ConfigWaveForm(min_period=PD["bounds"][0],
                                                max_period=PD["bounds"][1],
                                                taper_type="hann",
                                                taper_percentage=0.15)
        elif PD["adj_src_type"] == "cc_traveltime_misfit":
            cfgout = pyadjoint.ConfigCrossCorrelation(
                                                    min_period=PD["bounds"][0],
                                                    max_period=PD["bounds"][1],
                                                    taper_type='hann',
                                                    taper_percentage=0.3,
                                                    measure_type='dt',
                                                    use_cc_error=True,
                                                    dt_sigma_min=1.0,
                                                    dlna_sigma_min=0.5)
        elif PD["adj_src_type"] == "multitaper_misfit":
            cfgout = pyadjoint.ConfigMultiTaper(min_period=PD["bounds"][0],
                                                max_period=PD["bounds"][1],
                                                lnpt=15,
                                                transfunc_waterlevel=1e-10,
                                                water_threshold=0.02,
                                                ipower_costaper=10,
                                                min_cycle_in_window=0.5,
                                                taper_type='hann',
                                                taper_percentage=0.3,
                                                mt_nw=4.0,
                                                num_taper=5,
                                                dt_fac=2.0,
                                                phase_step=1.5,
                                                err_fac=2.5,
                                                dt_max_scale=3.5,
                                                measure_type='dt',
                                                dt_sigma_min=1.0,
                                                dlna_sigma_min=0.5,
                                                use_cc_error=True,
                                                use_mt_error=False)
        else:
            raise Exception("Pyadjoint 'adj_src_type' incorrectly specified")

    return cfgout
    
# ================================ RUN SCRIPTS =================================
def run_pyflex(PD,st,inv,event):
    """use pyflex to grab windows, current config set to defaults found on docs
    if writing into pyasdf files, window objects must be deconstructed and fed
    in as a parameter dictionary

    :type PD: dictionary
    :param PD: parameter sets
    :rtype windows: pyflex.window object
    :return windows: windows containing selected timespans where waveforms
    have acceptable match, also contains information about the window (see docs)
    """
    if PD["verbose"]:
        print("Running pyflex for [{}] configuration".format(
                                                        PD["pyflex_config"]))

    CD = choose_config("pyflex",PD)
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
    # create stalta data from envelopes of synthetic data
    windows,staltas = {},{}
    empties = 0
    for comp in PD["comp_list"]:
        print(comp,end='... ')

        obs,syn = breakout_stream(st.select(component=comp))
        window = pyflex.select_windows(observed=obs,
                                        synthetic=syn,
                                        config=config,
                                        event=pf_event,
                                        station=pf_station,
                                        plot=False)

        # calculate stalta
        syn_envelope = envelope(syn[0].data)
        stalta = pyflex.stalta.sta_lta(data=syn_envelope,
                                       dt=syn[0].stats.delta,
                                       min_period=PD["bounds"][0])
        staltas[comp] = stalta

        # check if pyflex is returning empty windows
        print("{} window(s)".format(len(window)))
        if not window:
            empties+=1
            continue
        windows[comp] = window

    # if all components show empty windows, raise the alarm
    if empties == len(PD["comp_list"]):
        raise Exception("Empty windows")

    # append STA/LTA water level to Par. dict. for plotting
    PD["stalta_wl"] = CD[0]

    # save windows into pyasdf file with stalta as the data and window-
    # parameter dictionaries as external information. dictionary needs
    # to be modified to work in pyasdf format
    if PD["dataset"]:
        if PD["verbose"]:print("Saving windows to PyASDF dataset")
        for comp in windows.keys():
            for i,window in enumerate(windows[comp]):
                internalpath = "{net}/{sta}_{comp}_{i}_{m}".format(
                                                            evid=PD["event_id"],
                                                            net=PD["network"],
                                                            sta=PD["station"],
                                                            m=PD["model"]
                                                            comp=comp,
                                                            i=i)
                # auxiliary data requires a data object, even though we only
                # want the window parameter dictionary. to save on space
                winnDixie = create_window_dictionary(window)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    PD["dataset"].add_auxiliary_data(data=np.array([True]),
                                                     data_type="MisfitWindows",
                                                     path=internalpath,
                                                     parameters=winnDixie)

    return windows, staltas, PD


def run_pyadjoint(st,windows,PD):
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
    if PD["verbose"]:
        print("Running pyAdjoint for type [{}] ".format(PD["adj_src_type"]))

    obs,syn = breakout_stream(st)
    net,sta,loc,cha = PD["code"].split('.')

    delta = st[0].stats.delta
    cfg = choose_config("pyadjoint",PD)

    # iterate through available window components
    adjoint_sources = {}
    for key in windows:

        obs_adj = obs.select(component=key)[0]
        syn_adj = syn.select(component=key)[0]

        # collect all windows into a single list object
        adjoint_windows = []
        for win in windows[key]:
            adj_win = [win.left*delta,win.right*delta]
            adjoint_windows.append(adj_win)

        adj_src = pyadjoint.calculate_adjoint_source(
                                             adj_src_type=PD["adj_src_type"],
                                             observed=obs_adj,
                                             synthetic=syn_adj,
                                             config=cfg,
                                             window=adjoint_windows,
                                             plot=False
                                             )
        adjoint_sources[key] = adj_src
                
        if PD["dataset"]:
            if PD["verbose"]:print("Saving adj src [{}] to pyASDF".format(key))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                adj_src.write_to_asdf(PD["dataset"],time_offset=0)
    
        if PD["save_adj_src"][0]:
            _check_path(join(PD["save_adj_src"][1],PD['event_id']))
            fidout = "{evid}/adjsrc_{net}_{sta}_{comp}_{mod}".format(
                                                        evid=PD['event_id'],
                                                        net=PD['network'],
                                                        sta=PD['station'],
                                                        mod=PD['model'],
                                                        comp=key)
            outpath = join(PD["save_adj_src"][1],fidout)
            adj_src.write(outpath,format="SPECFEM")

    return adjoint_sources


def bob_the_builder():
    """main processing script

    ++intra-function parameters and choices:
    :type EVENT_ID: list of str
    :param EVENT_ID: GEONET event ID's
    :type STANET_NAMES: list of str
    :param STANET_NAMES: identifiers for data fetching, in the form NN.SSS(S)
                         examples: NZ.HIZ, XX.RD01
    :type MINIMUM/MAXIMUM_FILTER_PERIOD: int
    :param MINIMUM/MAXIMUM_FILTER_PERIOD: bandpass filtering bounds in seconds
    :type ROTATE_TO_RTZ: bool
    :param ROTATE_TO_RTZ: RTZ if True, NEZ if False
    :type UNIT_OUTPUT: str
    :param UNIT_OUTPUT: from obspy, choices: "DISP", "VEL", "ACC"
    :type PYFLEX_CONFIG: str
    :param PYFLEX_CONFIG: choice of configuration for pyflex, choices:
                    "default": what is listed on the pyflex docs
                    "UAF": values taken from Tape et al. at Fairbanks
                    "NZ": values deemed best for New Zealand context
    :type ADJOINT_SRC_TYPE: str
    :param ADJOINT_SRC_TYPE: from pyadjoint, misfit calculation style, choices:
                    "cc_traveltime_misfit": squared traveltime differences
                    "multitaper_misfit": freq. dependent phgase differences
                    "waveform_misfit": squared difference of waveforms
    :type PLOT_*: bool
    :param PLOT_*: plot waveform (WAV), src-rcv map (MAP), active faults take 
                   some time to plot, faster without (FAULTS_ON_MAP)
    :type SAVE_*: tuple -> (bool,str)
    :param SAVE_*: save figures (PLOT) or all data (PYASDF), if so output path
                   adjoint sources can be saved on their own for easy input 
                   to specfem (ADJSRC_SEPARATE)
    :type VERBOSE: bool
    :param VERBOSE: enable print statements throughout 
    """
    
    ALLSTATIONS = {"GEONET":['NZ.BFZ','NZ.BKZ','NZ.HAZ','NZ.HIZ','NZ.KNZ',
                    'NZ.MRZ','NZ.MWZ','NZ.OPRZ','NZ.PUZ','NZ.PXZ','NZ.RTZ',
                    'NZ.TLZ','NZ.TOZ','NZ.TSZ','NZ.VRZ','NZ.WAZ'],
                   "FATHOM":['XX.RD01','XX.RD02','XX.RD03','XX.RD04','XX.RD05',
                    'XX.RD06','XX.RD07','XX.RD08','XX.RD09','XX.RD10',
                    'XX.RD11','XX.RD12','XX.RD13','XX.RD14','XX.RD15',
                    'XX.RD16','XX.RD17','XX.RD18','XX.RD19','XX.RD20',
                    'XX.RD21','XX.RD22'],
                    "TEST":['NZ.BKZ']
                    }
    # ============================ vPARAMETER SETv =============================
    # MODEL
    MODEL_NUMBER = 0
    # SOURCE-RECEIVER
    EVENT_IDS = ["2014p240655"]
    STANET_CHOICE = "FATHOM"
    # >> PREPROCESSING
    MINIMUM_FILTER_PERIOD = 6
    MAXIMUM_FILTER_PERIOD = 30
    ROTATE_TO_RTZ = True
    UNIT_OUTPUT = "DISP"
    # >> PYFLEX
    PYFLEX_CONFIG = "UAF"
    # >> PYADJOINT
    ADJOINT_SRC_TYPE = "multitaper_misfit"
    # >> PLOTTING
    PLOT_WAV = True
    PLOT_MAP = False
    PLOT_FAULTS_ON_MAP = True
    SHOW_PLOTS = True
    # >> SAVING
    SAVE_PLOT = (False,pathnames()["adjtomoplots"])
    SAVE_PYASDF = (True,pathnames()["adjtomodata"] + "PYASDF")
    SAVE_ADJSRC_SEPARATE = (False,pathnames()["adjtomodata"] + "ADJSRC")
    # >> MISC.
    VERBOSE = True
    LOG = (False,pathnames()["adjtomodata"] + "LOGS")
    # ============================ ^PARAMETER SET^ =============================
    # START LOGGING TO TEXTFILE
    if LOG[0]:
        sys.stdout = open(join(LOG[1],_create_log_fid(LOG[1])),"w")
    
    # PARAMETER AUTO SET
    COMPONENT_LIST = ["N","E","Z"]
    if ROTATE_TO_RTZ:
        COMPONENT_LIST = ["R","T","Z"]
    STANET_NAMES = ALLSTATIONS[STANET_CHOICE]    
    MODEL_NUMBER = "m{:0>2}".format(MODEL_NUMBER)
    
    # PRINT PARAMETERS
    if VERBOSE:
        # print the parameters in std. out
        import time
        template = ("\nPARAMETERS {time}\n{lines}\n\n"
                    "Model:       {mod}\n"
                    "Stations:    {sta}\n"
                    "Bandpass:    {tmin},{tmax}\n"
                    "Rotate:      {rot}\n"
                    "Units:       {uni}\n"
                    "Pyflex:      {pyf}\n"
                    "Pyadjoint:   {pya}\n"
                    "Plot:        Wav={wav}, Map={map}, Faults={fau}\n"
                    "Save plot:   {sav0} to {sav1}\n"
                    "Show plot:   {show}\n"
                    "Save Pyasdf: {asdf0} to {asdf1}\n"
                    "Verbose:     {ver}\n{lines}"
                    )
        print(template.format(mod=MODEL_NUMBER,sta=STANET_CHOICE,
                              tmin=MINIMUM_FILTER_PERIOD,
                              tmax=MAXIMUM_FILTER_PERIOD,rot=ROTATE_TO_RTZ,
                              uni=UNIT_OUTPUT,pyf=PYFLEX_CONFIG,
                              pya=ADJOINT_SRC_TYPE,wav=PLOT_WAV,map=PLOT_MAP,
                              fau=PLOT_FAULTS_ON_MAP,sav0=SAVE_PLOT[0],
                              sav1=SAVE_PLOT[1],show=SHOW_PLOTS,
                              asdf0=SAVE_PYASDF[0],asdf1=SAVE_PYASDF[1],
                              ver=VERBOSE,lines="_"*75,
                              time=UTCDateTime()))
        time.sleep(3)

    # MAIN ITERATE OVER EVENTS
    for EVENT_ID in EVENT_IDS:
        print("===={}====".format(EVENT_ID))

        # if data should be stored, initiate pyasdf dataset, also reads
        # existing pyasdf datasets if they already exist
        if SAVE_PYASDF[0]:
            datasetname = join(SAVE_PYASDF[1],EVENT_ID+'.h5')
            DATASET = pyasdf.ASDFDataSet(datasetname,compression="gzip-3")
        else:
            DATASET = None

        for STANET_NAME in STANET_NAMES:
            print("\n{}\n".format(STANET_NAME))
            PAR_DICT = {"model":MODEL_NUMBER,
                        "network":STANET_NAME.split('.')[0],
                        "station":STANET_NAME.split('.')[1],
                        "code":"{}.*.HH?".format(STANET_NAME),
                        "event_id":EVENT_ID,
                        "bounds":(MINIMUM_FILTER_PERIOD,
                                  MAXIMUM_FILTER_PERIOD),
                        "rotate":ROTATE_TO_RTZ,
                        "units":UNIT_OUTPUT,
                        "pyflex_config":PYFLEX_CONFIG,
                        "adj_src_type":ADJOINT_SRC_TYPE,
                        "save_adj_src":SAVE_ADJSRC_SEPARATE,
                        "comp_list":COMPONENT_LIST,
                        "save_plot":SAVE_PLOT,
                        "plot":(PLOT_WAV,PLOT_MAP,PLOT_FAULTS_ON_MAP),
                        "show":SHOW_PLOTS,
                        "dataset":DATASET,
                        "verbose":VERBOSE
                        }

            # MAIN PROCESSING
            try:
                st,inv,event = initial_data_gather(PAR_DICT)
                windows,staltas,PAR_DICT = run_pyflex(PAR_DICT,st,inv,event)
                adj_src = run_pyadjoint(st,windows,PAR_DICT)
                build_figure(st,inv,event,windows,staltas,adj_src,PAR_DICT)
            except KeyboardInterrupt:
                sys.exit('Keyboard Interrupt')
            except Exception as e:
                print(e)
                continue


# =================================== MAIN ====================================
if __name__ == "__main__":
    bob_the_builder()
    # func_test.test_build_figure()
