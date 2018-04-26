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
import pytz
import glob
import traceback
import numpy as np
from obspy import read, read_inventory, Stream
from obspy.signal.cross_correlation import correlate

# internal packages
sys.path.append("../modules")
from getdata import pathnames

from utils import check_save, create_min_max, already_processed
from plotutils import plot_arrays, stacked_plot

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


# ============================ PROCESSING HELPERS ==============================

def preprocess(st_raw,inv,resample,night=False,water_level=60):
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
    # trim traces to nz local night time, NZ ~= UTC+12
    if night:
        sunset_in_nz = st_pp[0].stats.starttime + 6*3600
        sunrise_in_nz = st_pp[0].stats.starttime + 18*3600
        st_pp.trim(sunset_in_nz,sunrise_in_nz)

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

def create_TEORRm_arrays(st_raw, inv, night, already_preprocessed=False):
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
    if not already_preprocessed:
        st = preprocess(st_raw,inv,resample=sampling_rate_Hz,night=night)
    else:
        st = st_raw
    try:
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


        # determine median values for amplitude ratio by 5min and by 1hour
        Rm,Rh = [],[]
        minute_,hour_ = sampling_rate_min*5,sampling_rate_min*60
        for i,time_window in enumerate([minute_,hour_]):
            for S0 in range(0,len(amplitude_ratio)-time_window,time_window):
                S1 = S0 + time_window
                avg = np.median(amplitude_ratio[S0:S1])
                Rm.append(avg) if i==0 else Rh.append(avg)

        print(round(time.time()-T0,2))

        TEORRm = {"T":np.array(tremor_horizontal),
                  "E":np.array(earthquake_horizontal),
                  "O":np.array(ocean_horizontal),
                  "R":np.array(amplitude_ratio),
                  "Rm":np.array(Rm),
                  "Rh":np.array(Rh)
                  }

        return st, TEORRm

    except Exception as e:
        print('[TEORRm creation error]')
        return st, None

def data_gather_and_process(code_set,pre_filt=False,night=False):
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
    net,sta,loc,cha,d,year,jday = code_set.split('.')

    fid_path = pathnames()['RDF'] + "{y}/XX/{s}/HH{c}.D".format(y=year,
                                                              s=sta,
                                                              c="{c}")
    inv_path = pathnames()['RDF'] + "XX.RDF.DATALESS"

    # check what combination of files are available
    path_dict = check_save(code_set,night=night)
    process_check = False

    if not path_dict["npz"]:
        if not path_dict["pickle"]:
            print("[files don't exist, processing...]")
            st = Stream()
            inv = read_inventory(inv_path)
            for comp in ["N","E"]:
                fid = os.path.join(fid_path,code_set).format(c=comp)
                st += read(fid)
            if len(st) != 2:
                print("{} traces found, skipping".format(len(st)))
                return False,False,False
        else:
            inv = None
            process_check = True
            st = read(path_dict['pickle'])

        st_proc, TEORRm = create_TEORRm_arrays(st,inv,night=night,
                                        already_preprocessed=process_check)
        check_bool = check_save(code_set,st=st_proc,TEORRm=TEORRm,night=night)
        if not check_bool:
            return False,False,False

    else:
        print("[files exist, reading...]")
        st_proc = read(path_dict['pickle'])
        TEORRm = np.load(path_dict['npz'])

    if pre_filt:
        st_proc.filter('bandpass',freqmin=pre_filt[0],freqmax=pre_filt[1])

    # count tremors
    sig_arrays = tremor_counter(TEORRm["Rm"],avg_choice="median",night=night)

    return st_proc, TEORRm, sig_arrays

def __mean_std_creator(night=False):
    """determine mean value of traces for all TEORRm arrays, and one-sigma value
    to be used for counting tremors. nighttime only looks at files with suffix
    _night
    """
    filepath = pathnames()['data'] + 'TEROR/'
    N = "_night"
    if night:
        N = ""
    terorfiles = glob.glob(filepath + 'XX*{N}.npz'.format(N=N))

    # file check
    num_of_files = len(terorfiles)
    outpath = os.path.join(filepath,'Rm_mean_std{N}_{C}.npz'.format(N=N,
                                                                C=num_of_files))
    if os.path.exists(outpath):
        print("[mean/std file exists, reading...]")
        avgs_and_one_sigma = np.load(outpath)
        mean_ = avgs_and_one_sigma['mean']
        median_ = avgs_and_one_sigma['median']
        one_sigma_ = avgs_and_one_sigma['onesigma']

        return mean_,median_,one_sigma_

    # if new file, calculate all averages etc.
    means_,medians_,one_sigmas_ = [],[],[]
    for count,fid in enumerate(terorfiles):
        TEORRm = np.load(fid)
        Rm = TEORRm['Rm']
        means_.append(np.mean(Rm))
        medians_.append(np.median(Rm))
        one_sigmas_.append(np.std(Rm))

    mean_ = np.mean(means_)
    median_ = np.mean(medians_)
    one_sigma_ = np.mean(one_sigmas_)

    print("[creating mean/std file...]\n\t"
            "Mean:{mean},Med:{med},1-sig:{one_sig}".format(
                                                    mean=round(mean_,2),
                                                    med=round(median_,2),
                                                    one_sig=round(one_sigma_,2)
                                                    ))

    np.savez(outpath,mean=mean_,median=median_,onesigma=one_sigma_)

    # delete any previous files that are superceded by current write
    outpath_wild = os.path.join(filepath,'Rm_mean_std{N}_{C}.npz'.format(N=N,
                                                                        C="*"))
    files_to_delete = glob.glob(outpath_wild)
    files_to_delete.remove(outpath)
    if files_to_delete:
        for fid in files_to_delete:
            os.remove(fid)

    return mean_,median_,one_sigma_

def tremor_counter(Rm,avg_choice='mean',night=False):
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
    mean_,median_,one_sigma = __mean_std_creator(night)

    avg_val = mean_
    if avg_choice == 'median':
        avg_val = median_

    two_sigma = avg_val + one_sigma * 2
    three_sigma = avg_val + one_sigma * 3

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
    return Rm_2sig, Rm_3sig, two_sigma


def convert_UTC_to_local(st,local_timezone="Pacific/Auckland"):
    """convert timezone of a stream object from UTC to local timezone
    """
    # create UTC aware datetime objects
    startUTC = st[0].stats.starttime.datetime.replace(tzinfo=pytz.utc)
    endUTC = st[0].stats.endtime.datetime.replace(tzinfo=pytz.utc)
    startLOC = startUTC.astimezone(pytz.timezone(local_timezone))
    endLOC = endUTC.astimezone(pytz.timezone(local_timezone))

    return startLOC, endLOC

def __handler(signum,frame):
    raise Exception("[process timed out]")


# ============================= MAIN PROCESSING ================================
def stacked_process(jday):
    """creating stacked plots by running processing for multiple stations and
    feeding the outputs into plotting script
    data manipualation includes:
    -creating min/max arrays to ease plotting requirements
    -normalizing traces -1 to 1 to allow plotting on the same figure
    -removing values greater than 0.5 (most likely earthquake signals)
    """
    # ///////////////////// parameter set \\\\\\\\\\\\\\\\\\\\\\\
    station_list = [8,9,12,13,14,16,6,1]
    nighttime_only = True
    stop_if_tremor_num_below = 3
    stop_if_stations_above = len(station_list) // 2
    # \\\\\\\\\\\\\\\\\\\\\ parameter set ///////////////////////

    # accumulate, process all data and place in arrays for later plotting
    code_set_template = "XX.RD{s:0>2}.10.HH{c}.D.2017.{d}"
    num_low_tremor_events=0

    y_N_list,y_E_list,Rm_list,sig_list,tremor_list = [],[],[],[],[]
    ano_NE_list = [[],[]]

    for station in station_list:
        code_set = code_set_template.format(s=station,c="{c}",d=jday)
        print('\n',code_set)
        try:
            st_,TEORRm,sig_arrays = data_gather_and_process(
                                code_set,pre_filt=[2,8],night=nighttime_only)
            if not st_:
                continue
        except Exception as e:
            print("ERROR {jday} RD{station}".format(jday=jday,station=station))
            traceback.print_exc()
            continue

        # hacky way to avoid errors with overwriting the last entry before plot
        st = st_
        sig2_array,sig3_array,two_sigma = sig_arrays

        # check stop if detection threshold met
        if stop_if_tremor_num_below:
            if len(sig2_array[sig2_array>0]) < stop_if_tremor_num_below:
                print("Low tremor event {0}/{1}".format(num_low_tremor_events+1,
                                                      stop_if_tremor_num_below))
                num_low_tremor_events += 1
                if (num_low_tremor_events > stop_if_stations_above):
                    print("Exiting due to low detection")
                    return False
                continue

        # set up waveform plotting arrays
        for comp in ["N","E"]:
            x,y = create_min_max(st.select(component=comp)[0])
            y_E_max = y.max()
            if comp == "N":
                y_N_max = y.max()
            y/=y.max()
            y_N_list.append(y) if comp == "N" else y_E_list.append(y)

        # create annotations from max amplitude values
        for i,max_ in enumerate([y_N_max,y_E_max]):
            ano_NE_list[i].append("{s} {a}um/s".format(s=code_set.split('.')[1],
                                                        a=round(max_ * 1E6,2)
                                                        ))

        # normalize sigma and Rm 0to1
        ratio_median = TEORRm["Rm"]
        # ratio_median/=ratio_median.max()
        Rm_list.append(ratio_median)

        # if not sig2_array.max() == 0:
        #     sig2_array/=sig2_array.max()
        sig2_array[sig2_array==0] = np.nan
        sig_list.append(sig2_array)

    if not y_N_list:
        return False

    # RMS of tremor signal for envelope plots, nan out amplitudes above 0.5
    for N,E in zip(y_N_list,y_E_list):
        horizontal_rms = 0.5 * (N**2 + E**2)
        horizontal_rms /= horizontal_rms.max()
        for data in N,E,horizontal_rms:
            data[abs(data)>0.5] = np.nan

        tremor_list.append(horizontal_rms)


    # time arrays
    startNZ,endNZ = convert_UTC_to_local(st)
    t0 = startNZ.hour
    t1 = t0 + 12
    t = np.linspace(t0,t1,len(x))

    stacked_plot(code=code_set,x=t,N=y_N_list,E=y_E_list,
                 Rm_list=Rm_list,
                 sig_list=sig_list,
                 ano_list=ano_NE_list,
                 tremor_list=tremor_list,
                 night=nighttime_only,
                 horizontal_line=two_sigma,
                 show=True,
                 save=False)

    return True

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
        plot_arrays(st,code_set,TEORRm,sig,show=True,save=True)


if __name__ == "__main__":
    for jday in range(262,300):
        print('=========={}========='.format(jday))
        stacked_process(jday)
