"""Python port of the modified Frequency Scanning Method (mFSM),
from Katakami et al. 2017 JGR
Used to detect tremors in waveform data by taking amplitude ratio of different
frequency filter bands, removes earthquake data using a simple cross-correlation
with an exponential function. Counts tremor detection using standard deviation
thresholds.
"""
import os
import sys
import time
import numpy as np
from obspy import read, read_inventory, Stream
from obspy.signal.cross_correlation import correlate

# internal packages
sys.path.append("../modules")
from getdata import pathnames

from utils import z2nan, check_save, create_min_max
from plotutils import plot_arrays, stacked_plot

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


def preprocess(st_raw,inv,resample,water_level=60):
    """preprocess waveform data: resample, detrend, taper, remv. resp.
    :type st_raw: obspy stream
    :param st_ray: raw data stream
    :type inv: obspy inventory
    :param inv: response information for station
    :type resample: int
    :param resample: new sampling rate to set
    :type water_level: int
    :param water_level: water level input during instrument response
    """
    print("[preprocess]",end=" ")
    T0 = time.time()
    st_pp = st_raw.copy()
    st_pp.resample(resample)
    st_pp.detrend("demean") # not in original code
    st_pp.detrend("linear") # not in original code
    st_pp.taper(max_percentage=0.05)
    st_pp.attach_response(inv)
    pre_filt = [.001,.01,50,55]
    st_pp.remove_response(output="VEL",
                          pre_filt=pre_filt, # not in original code
                          water_level=water_level, # not in original code
                          plot=False)
    print(round(time.time()-T0,2))
    return st_pp

def create_horizontal_data(st,bounds):
    """filter north and east component data and combine into horizontal comp
    :type st: obspy stream
    :param st: preprocessed data
    :type bounds: list of floats
    :param bounds: [lower bound,upper bound] for filtering
    :rtype data_horizontal: numpy array
    :return data_horizontal: average of horizontal components
    """
    st_chd = st.copy()
    st_chd.filter("bandpass",freqmin=bounds[0],
                             freqmax=bounds[1],
                             corners=3,
                             zerophase=True)
    data_north = st_chd.select(component="N")[0].data
    data_east = st_chd.select(component="E")[0].data

    data_horizontal = np.sqrt(data_north**2 + data_east**2)

    return data_horizontal

def set_water_level(st,band):
    """create water level from ocean band filter bands to avoid division by
    very small values, which could lead to false detections
    :type st: obspy stream
    :param st: preprocessed data
    :type band: list of floats
    :param band: [lower bound,upper bound] for filtering
    :rtype data_horizontal: numpy array
    :return data_horizontal: average of horizontal components, with WL set
    :rtype water_level: float
    :return water_level: mean value of data_horizontal
    """
    data_horizontal = create_horizontal_data(st,band)
    water_level = data_horizontal.mean()
    low_values_flags = data_horizontal < water_level
    data_horizontal[low_values_flags] = water_level

    return data_horizontal, water_level

def detect_earthquakes(tremor_horizontal,sampling_rate_min,corr_criteria=0.7):
    """remove regular earthquake from waveforms by taking correlations with
    an exponential function. If correlation criteria met, earthquake 'detected'
    :type tremor_horizontal: numpy array
    :param tremor_horizontal: average of horizontal components in Tremor band
    :type sampling_rate_min: float
    :param sampling_rate_min: number of samples in one minute
    :type corr_criteria: float
    :param corr_criteria: threshold for detecting earthquakes, defaults to 0.7
    :rtype quakearray: numpy array
    :return quakearray: array containing -1's for detected earthquakes
    """
    print("[detect_earthquakes]",end=" ")
    T0 = time.time()
    x= np.linspace(0.002,6,sampling_rate_min)
    exp_internal = -(x/2)*2
    exp_template = np.exp(exp_internal)

    # set different number of samples for array lenghts: .5-minute, 1.5-minute
    sampling_rate_half_min = int(sampling_rate_min * (1/2))
    sampling_rate_one_one_half_min = int(sampling_rate_min * (3/2))

    # scan through datastream by time windows of length sampling_rate_min,
    # use a full window shifting cross correlation
    quakearray = np.array([])

    # fill-value arrays if earthquake detected
    neg_ones = -1*(np.ones(sampling_rate_min))
    neg_ones_ext = -1*(np.ones(sampling_rate_one_one_half_min))

    endtrace = len(tremor_horizontal)
    quakecount = 0
    for S0 in range(0,endtrace,sampling_rate_min):
        # print("{0}/{1}".format(S0/sampling_rate_min,endtrace/sampling_rate_min))
        S1 = S0 + sampling_rate_min

        tremor_snippet = tremor_horizontal[S0:S1]
        exp_correlation = correlate(exp_template,tremor_snippet,shift=len(x))
        max_corr = exp_correlation.max()
        # if correlation criteria met
        if max_corr > corr_criteria:
            quakecount +=1
            if S0 == 0:
                quakearray = np.append(quakearray,neg_ones)
            else:
                quakearray_new = quakearray[:S0-sampling_rate_half_min]
                quakearray = np.append(quakearray_new,neg_ones_ext)
        else:
            quakearray = np.append(quakearray,tremor_snippet)

    print(round(time.time()-T0,2))
    print("{} earthquakes detected".format(quakecount))

    return quakearray

def create_TEORRm_arrays(st_raw, inv):
    """create filtered horizontal bands, set water level, perform simple
    earthquake detection and create amplitude ratios.
    :type st_raw: obspy stream
    :param st_raw: raw stream
    :type inv: obspy inventory
    :param inv: response information for st_raw
    :rtpye TEORRm: list of numpy arrays
    :return TEORRm: arrays containing filtered waveforms and amp. ratios
    """
    # parameter set
    tremor_band = [2,8]
    earthquake_band = [10,20]
    ocean_band = [0.5,1]

    sampling_rate_Hz = 50
    sampling_rate_min = 50*60

    # preprocess
    st = preprocess(st_raw,inv,resample=sampling_rate_Hz)

    # create arrays for different freq bands
    tremor_horizontal = create_horizontal_data(st,tremor_band)
    earthquake_horizontal = create_horizontal_data(st,earthquake_band)

    # set water level on ocean band
    ocean_horizontal, water_level = set_water_level(st,ocean_band)

    # simple earthquake detection
    quakearray = detect_earthquakes(tremor_horizontal,sampling_rate_min)
    quakearray_mean = quakearray.mean()

    amplitude_ratio = []
    ratio_equation = lambda T,E,O: T**2 / (E*O)
    print("[amplitude_ratios]",end=" ")
    T0 = time.time()
    for S0,_ in enumerate(tremor_horizontal):
        if quakearray[S0] == -1:
            amplitude_ratio.append(-1)
        else:
            R = ratio_equation(T=tremor_horizontal[S0],
                               E=earthquake_horizontal[S0],
                               O=ocean_horizontal[S0]
                               )
            amplitude_ratio.append(R)

    print(round(time.time()-T0,2))

    # determine median values for amplitude ratio
    print("[median values]",end=" ")
    time_window = sampling_rate_min * 5
    median_amp_ratio = []
    for S0 in range(0,len(amplitude_ratio)-time_window,time_window):
        S1 = S0 + time_window
        med = np.median(amplitude_ratio[S0:S1])
        median_amp_ratio.append(med)

    print(round(time.time()-T0,2))

    TEORRm = [np.array(tremor_horizontal),
              np.array(earthquake_horizontal),
              np.array(ocean_horizontal),
              np.array(amplitude_ratio),
              np.array(median_amp_ratio)
              ]

    return st, TEORRm

def data_gather_and_process(code_set,pre_filt=False):
    """grab relevant data files for instrument code, process using internal
    functions, return arrays containing filtered waveforms and ratio values
    :type code_set: str
    :param code_set: instrument code set in main
    :rtype TEORRm: list of numpy arrays
    :return TEORRm: arrays containing filtered waveforms and amp. ratios
    :rtype sig: list of numpy arrays
    :return sig: two column list, 2-sigma and 3-sigma of Rm respectively
    """
    # setting up datapaths
    net,sta,loc,cha,year,jday = code_set.split('.')

    fid_path = pathnames()['RDF'] + "{y}/XX/{s}/HH{c}".format(y=year,
                                                              s=sta,
                                                              c="{c}")
    inv_path = pathnames()['RDF'] + "XX.RDF.DATALESS"

    path_dict = check_save(code_set)


    if not path_dict:
        print("[files don't exist, processing...]")
        # read in data
        st = Stream()
        inv = read_inventory(inv_path)
        for comp in ["N","E"]:
            fid = os.path.join(fid_path,code_set).format(c=comp)
            st += read(fid)
        # process
        st_proc, TEORRm = create_TEORRm_arrays(st,inv)
        _ = check_save(code_set,st=st_proc,TEORRm=TEORRm)
    else:
        print("[files exist, reading...]")
        st_proc = read(path_dict['pickle'])
        TEORRm = np.load(path_dict['npz'])
        TEORRm = [TEORRm['T'],TEORRm['E'],TEORRm['O'],TEORRm['R'],TEORRm['Rm']]

    # filter streams for plotting
    if pre_filt:
        st_proc.filter('bandpass',freqmin=pre_filt[0],freqmax=pre_filt[1])

    # count tremors
    Rm = TEORRm[-1]
    Rm_2sig,Rm_3sig = tremor_counter(Rm)
    sig = [Rm_2sig,Rm_3sig]

    return st_proc, TEORRm, sig

def tremor_counter(Rm,nighttime=False):
    """port of calc_numTT from Satoshi, counts the number of tremors per day
    using standard deviations to determine tremor threshold, for 2-sigma
    and 3-sigma detection thresholds
    :type Rm: numpy array
    :param Rm: array containing the median values for amplitude ratios
    :type nighttime: boolean (not implemented yet)
    :param nighttime: investigate only nighttime signals to avoid cultural noise
    :rtype Rm_?sig: numpy array
    :return Rm_?sig: ?-sigma detection array for Rm
    """
    one_sigma = np.std(Rm[Rm>0])
    mean_val = np.mean(Rm[Rm>0])
    two_sigma = mean_val + one_sigma * 2
    three_sigma = mean_val + one_sigma * 3

    Rm_2sig = np.copy(Rm)
    Rm_3sig = np.copy(Rm)

    # count tremor activity by fulfilling threshold criteria for 2sig and 3sig
    for i,section in enumerate(Rm):
        if section < two_sigma:
            Rm_2sig[i] = 0
            Rm_3sig[i] = 0
        elif section < three_sigma:
            Rm_3sig[i] = 0

    print("{} tremors detected at 3-sigma detection".format(
                                                    len(Rm_3sig[Rm_3sig>0]))
                                                    )
    print("{} tremors detected at 2-sigma detection".format(
                                                    len(Rm_2sig[Rm_2sig>0]))
                                                    )
    return Rm_2sig, Rm_3sig

def time_convert():
    """convert time from UTC to local-NZ time. for use in e.g. finding out night
    time to remove the effect of cultural noise
    """

# ============================= MAIN PROCESSING ================================
def stacked_process():
    """creating stacked plots by running processing for multiple stations and
    feeding the outputs into plotting script
    data manipualation includes:
    -creating min/max arrays to ease plotting requirements
    -normalizing traces -1 to 1 to allow plotting on the same figure
    -removing values greater than 0.5 (most likely earthquake signals)
    """
    # ///////////////////// parameter set \\\\\\\\\\\\\\\\\\\\\\\
    station_list = [8,9,12,13,14]
    jday = 220
    # \\\\\\\\\\\\\\\\\\\\\ parameter set ///////////////////////

    # accumulate all data
    code_set_template = "XX.RD{s:0>2}.10.HH{c}.2017.{d}"
    y_N_list,y_E_list,Rm_list,sig_list,sta_list = [],[],[],[],[]
    for station in station_list:
        code_set = code_set_template.format(s=station,c="{c}",d=jday)
        sta_list.append(code_set.split('.')[1])
        st,TEORRm,sig = data_gather_and_process(code_set,pre_filt=[2,8])

        # set up plotting arrays
        for comp in ["N","E"]:
            x,y = create_min_max(st.select(component=comp)[0])
            y/=y.max()
            y[y>0.5]=np.nan
            if comp == "N":
                y_N_list.append(y)
            elif comp == "E":
                y_E_list.append(y)

        # grab sigma values and Rm values
        Rm = TEORRm[-1]
        Rm/=Rm.max()
        Rm_list.append(TEORRm[-1])

        sigma2 = sig[1]
        sigma2/=sigma2.max()
        sig_list.append(z2nan(sigma2))

    stacked_plot(x,y_N_list,y_E_list,Rm_list,sig_list,sta_list)

def single_process():
    """process a single station day by day
    """
    # ///////////////////// parameter set \\\\\\\\\\\\\\\\\\\\\\\
    code_set_template = "XX.RD13.10.HH{c}.2017.{d}"
    # \\\\\\\\\\\\\\\\\\\\\ parameter set ///////////////////////
    for i in range(220,250):
        code_set = code_set_template.format(c="{c}",d=i)
        print(code_set)
        st,TEORRm,sig = data_gather_and_process(code_set,pre_filt=[2,8])
        plot_arrays(st,code_set,TEORRm,sig,show=False,save=True)



if __name__ == "__main__":
    stacked_process()
