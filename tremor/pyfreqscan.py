"""Python port of the modified Frequency Scanning Method (mFSM),
from Katakami et al. 2017 JGR, with functionalities of FSM from Sit et al. 2012
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
from obspy import read, read_inventory, Stream, UTCDateTime
from obspy.signal.cross_correlation import correlate

# internal packages
sys.path.append("../modules")
from getdata import pathnames
from utils import check_save, create_min_max, already_processed, check_save
from plotutils import plot_arrays, stacked_plot, gridspec_plot

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# ========================== GENERAL HELPER FUNCTIONS ==========================

def convert_UTC_to_local(st,local_timezone="Pacific/Auckland"):
    """convert timezone of a stream object from UTC to local timezone
    :type st: obspy stream
    :param st: datastream with timezone assumed UTC
    :type local_timezone: string
    :param local_timezone: datetime timezone identifier
    :rtype *LOC: datetime.datetime
    :return *LOC: converted datetime objects
    """
    # create UTC aware datetime objects
    startUTC = st[0].stats.starttime.datetime.replace(tzinfo=pytz.utc)
    endUTC = st[0].stats.endtime.datetime.replace(tzinfo=pytz.utc)
    startLOC = startUTC.astimezone(pytz.timezone(local_timezone))
    endLOC = endUTC.astimezone(pytz.timezone(local_timezone))

    return startLOC, endLOC

def average_array(data,number_of_samples):
    """return a shortened version of input array which contains averaged
    values for every number_of_samples specified. the end member may contain a
    disproportionate number of samples to make up for the indivisibility of
    the number of samples by the chosen division rate
    :type data: np.array
    :param data: array to be averaged
    :type number_of_samples: int
    :param number_of_samples: length of window in samples to average over
    :rtype new_data: np.array
    :return new_data: averaged array
    """
    new_data = []
    for i in range(0,len(data)-number_of_samples,number_of_samples):
        j = i + number_of_samples
        if j > len(data)-number_of_samples:
            j = len(data) - 1
        new_data.append(np.mean(data[i:j]))

    return np.array(new_data)

def create_filt_horiz_data(st,bounds):
    """filter north and east component data and combine into horizontal comp,
    filter to given bounds
    :type st: obspy stream
    :param st: preprocessed data
    :type bounds: list of floats
    :param bounds: [lower bound,upper bound] for filtering
    :type take_average: int
    :param take_average: if averages of arrays should be taken, if so give an
    integer number of samples to take averaging windows over
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
    """create water level from filter bands to avoid division by small values
    :type st: obspy stream
    :param st: preprocessed data
    :type band: list of two floats
    :param band: [lower bound,upper bound] for filtering
    :rtype data_horizontal: np.array
    :return data_horizontal: average of horizontal components, with WL set
    :rtype water_level: float
    :return water_level: mean value of data_horizontal
    """
    data_horizontal = create_filt_horiz_data(st,band)
    water_level = data_horizontal.mean()
    low_values_flags = data_horizontal < water_level
    data_horizontal[low_values_flags] = water_level

    return data_horizontal, water_level

def mean_std_creator(night=False):
    """determine mean and 1-sigma values for all TEORRm arrays, to be used for
    counting tremors.
    :type night: bool
    :param night: consider only nighttime files with suffix '_night'
    :rtype sig_dict: dict
    :return sig_dict: contains mean, median and 1-sigma values for each ratio
    array Rs,Rm and Rh
    """
    filepath = os.path.join(pathnames()['data'],
                                        'TEROR','201?','*','*','npz','')
    outpath_template = os.path.join(pathnames()['data'],'TEROR',
                                                        'MEANSTD{N}_{C}.npz')
    N = "_night"
    if not night:
        N = ""
    terorfiles = glob.glob(filepath + 'XX*{N}.npz'.format(N=N))
    outpath = outpath_template.format(N=N,C=len(terorfiles))
    if os.path.exists(outpath):
        # print("[mean/std file exists, reading...]")
        sig_dict = np.load(outpath)
        return sig_dict

    # if new file, calculate all averages etc.
    sig_dict = {}
    for choice in ['Rs','Rm','Rh']:
        means_,medians_,one_sigmas_ = [],[],[]
        for count,fid in enumerate(terorfiles):
            TEORRm = np.load(fid)
            ratio = TEORRm[choice]
            means_.append(np.mean(ratio))
            medians_.append(np.median(ratio))
            one_sigmas_.append(np.std(ratio))

        # dictionary
        sig_dict["{}_mean".format(choice)] = np.mean(means_)
        sig_dict["{}_median".format(choice)] = np.mean(medians_)
        sig_dict["{}_sigma".format(choice)] = np.mean(one_sigmas_)

    print("[creating mean/std file...]")
    np.savez(outpath,**sig_dict)

    # delete any previous files that are superceded by current write
    outpath_wild = outpath_template.format(c=choice,N=N,C="*")
    files_to_delete = glob.glob(outpath_wild)
    files_to_delete.remove(outpath)
    if files_to_delete:
        for fid in files_to_delete:
            os.remove(fid)

    return sig_dict

# =========================== PROCESSING FUNCTIONS =============================

def preprocess(st_original,inv,resample,night=False,water_level=60):
    """preprocess waveform data: resample, detrend, taper, remv. resp.
    :type st_original: obspy stream
    :param st_original: raw data stream
    :type inv: obspy inventory
    :param inv: response information for station
    :type resample: int
    :param resample: new sampling rate to set
    :type water_level: int
    :param water_level: water level input during instrument response
    :rtype st_pp: obspy stream
    :return st_pp: preprocessed stream object
    """
    print("[preprocess]",end=" ")
    st_pp = st_original.copy()

    T0 = time.time()
    # trim traces to nz local night time, NZ ~= UTC+12
    if night:
        T = st_pp[0].stats.starttime
        sunset_in_nz = UTCDateTime(T.year,T.month,T.day,6,0,0)
        sunrise_in_nz = UTCDateTime(T.year,T.month,T.day,18,0,0)
        st_pp.trim(sunset_in_nz,sunrise_in_nz)

    # data quality checks
    if not st_pp:
        print("\t--streams fall outside trim")
        return False
    elif st_pp[0].stats.npts != st_pp[1].stats.npts:
        print("\t--streams have inconsistent npts {0} {1}".format(
                                                        st_pp[0].stats.npts,
                                                        st_pp[1].stats.npts))
        return False

    # standard preprocess steps
    decimate_by = int(st_pp[0].stats.sampling_rate//resample)
    st_pp.decimate(decimate_by)
    st_pp.detrend("demean")
    st_pp.detrend("linear")
    st_pp.taper(max_percentage=0.05)
    st_pp.attach_response(inv)
    st_pp.remove_response(output="VEL",
                          pre_filt=[.001,.01,50,55], # not in original code
                          water_level=water_level, # not in original code
                          plot=False)

    print(round(time.time()-T0,2),'s')

    return st_pp

def detect_earthquakes(de_array,sampling_rate,corr_criteria=0.7):
    """try to remove earthquake from waveforms by taking correlations with
    an exponential function. If correlation criteria met, earthquake 'detected'
    :type de_array: numpy array
    :param de_array: datastream representing waveform envelope
    :type sampling_rate: float
    :param sampling_rate: sampling rate
    :type corr_criteria: float
    :param corr_criteria: threshold for detecting earthquakes, defaults to 0.7
    :rtype quakearray: np.array
    :return quakearray: de_array containing -1's for detected earthquakes
    """
    T0 = time.time()

    # set exponential template
    sampling_rate_min = sampling_rate * 60
    sampling_rate_half_min = int(sampling_rate_min * (1/2))
    sampling_rate_one_one_half_min = int(sampling_rate_min * (3/2))
    x= np.linspace(0.002,6,sampling_rate_min)
    exp_internal = -(x/2)*2
    exp_template = np.exp(exp_internal)

    # fill-value arrays if earthquake detected
    nan_fill = np.nan*(np.ones(sampling_rate_min))
    nan_fill_ext = np.nan*(np.ones(sampling_rate_one_one_half_min))

    quakecount = 0
    quakearray = np.array([])
    for S0 in range(0,len(de_array),sampling_rate_min):
        S1 = S0 + sampling_rate_min
        tremor_snippet = de_array[S0:S1]
        exp_correlation = correlate(a=exp_template,
                                    b=tremor_snippet,
                                    shift=len(x))

        if exp_correlation.max() > corr_criteria:
            quakecount +=1
            if S0 == 0:
                quakearray = np.append(quakearray,nan_fill)
            else:
                quakearray_new = quakearray[:S0-sampling_rate_half_min]
                quakearray = np.append(quakearray_new,nan_fill_ext)
        else:
            quakearray = np.append(quakearray,tremor_snippet)

    print("[detect_earthquakes] {} quakes".format(quakecount),end=" ")
    print(round(time.time()-T0,2),'s')

    return quakearray

def tremor_counter(TEORRm,avg_choice="mean",stop_if_tremor_num_below=3,
                                                                night=False):
    """count number of tremors using a global 2-sigma threshold
    :type TEORRm: dict of np.arrays
    :param TEORRm: array containing filter-bands and ratios
    :type avg_choice: str
    :param avg_choice: mean or median for calculating 2 and 3-sigma
    :type stop_if_tremor_num_below: int
    :param stop_if_tremor_num_below: if number of tremors detected in Rm below,
    declare 'low tremor event' and skip
    :type night: bool
    :param night: consider only files with suffix '_night'
    """
    tremor_count_dict = {}
    for choice in ['Rs','Rm','Rh']:
        # set average values
        sig_dict = mean_std_creator(night)
        mean_ = sig_dict["{}_mean".format(choice)]
        median_ = sig_dict["{}_median".format(choice)]
        one_sigma = sig_dict["{}_sigma".format(choice)]

        avg_val = mean_
        if avg_choice == 'median':
            avg_val = median_

        two_sigma = avg_val + (2*one_sigma)
        three_sigma = avg_val + (3*one_sigma)

        # count tremors - can probably write cleaner
        R_2sig = np.copy(TEORRm[choice])
        R_3sig = np.copy(TEORRm[choice])
        for i,section in enumerate(TEORRm[choice]):
            if section < two_sigma:
                R_2sig[i] = 0
                R_3sig[i] = 0
            elif section < three_sigma:
                R_3sig[i] = 0

        print("[tremor_counter {c}: 2-sig:{two} / 3-sig:{three}]".format(
                                                    c=choice,
                                                    two=len(R_2sig[R_2sig>0]),
                                                    three=len(R_3sig[R_3sig>0])
                                                    ))


        # check-stop if low-tremor detection threshold met on Rs
        if stop_if_tremor_num_below and choice=='Rm':
            if len(R_2sig[R_2sig>0]) < stop_if_tremor_num_below:
                print("\t--tremor number below {n} ".format(
                                                n=stop_if_tremor_num_below))
                return None

        # nan out 0's for plotting
        R_2sig[R_2sig==0] = np.nan
        tremor_count_dict[choice] = R_2sig

    return tremor_count_dict

def create_modified_TEROR(st_raw, inv, night, already_preprocessed=False):
    """create passband arrays more in line with methods by Sit et al. 2012
    :type st_raw: obspy stream
    :param st_raw: raw stream
    :type inv: obspy inventory
    :param inv: response information for st_raw
    :rtpye TEORRm: list of numpy arrays
    :return TEORRm: arrays containing filtered waveforms and amp. ratios
    """
    # parameter set
    tremor_band = [2,8]
    tremor_band_short = [2,5]
    earthquake_band = [10,20]
    earthquake_band_short = [10,15]
    ocean_band = [0.5,1]
    surface_wave_band = [0.02,0.1]

    sampling_rate = 50
    five_seconds = sampling_rate * 5
    one_minute = sampling_rate * 60
    five_minutes = one_minute * 5
    one_hour = one_minute * 60

    # preprocess
    if not already_preprocessed:
        st = preprocess(st_raw,inv,resample=sampling_rate,night=night)
        if not st:
            return None, None
    else:
        st = st_raw
    try:
        # create arrays for different freq bands
        tremor_horizontal = create_filt_horiz_data(st,tremor_band_short)
        earthquake_horizontal = create_filt_horiz_data(st,earthquake_band)
        surface_wave_horizontal = create_filt_horiz_data(st,surface_wave_band)

        # simple earthquake detection and removal
        tremor_horizontal = detect_earthquakes(tremor_horizontal,sampling_rate,
                                                            corr_criteria=0.65)

        # create ratio arrays by time window, honor earthquake detection
        ratio_dict = {}
        ratio_equation = lambda T,E,O: T**2 / (E*O)
        for window,Rx in zip([five_seconds,five_minutes,one_hour],
                                                            ['Rs','Rm','Rh']):
            T = average_array(tremor_horizontal,window)
            E = average_array(earthquake_horizontal,window)
            S = average_array(surface_wave_horizontal,window)
            R = ratio_equation(T=T,E=E,O=S)
            R[np.isnan(R)] = 0
            ratio_dict[Rx] = R

        # 0 out any values of Rm > (Rs.mean + 3*Rs.1-sigma)
        Rh_mean_3sigma = ratio_dict['Rh'] + 3*np.std(ratio_dict['Rh'])
        Rh_m3s_long = np.repeat(Rh_mean_3sigma,
                                    len(ratio_dict['Rm'])//len(Rh_mean_3sigma))
        ratio_dict['Rm'][ratio_dict['Rm']>Rh_m3s_long] = 0

        # distribution
        TEROR = {"T":tremor_horizontal,
                  "E":earthquake_horizontal,
                  "S":surface_wave_horizontal,
                  "Rs":ratio_dict['Rs'],
                  "Rm":ratio_dict['Rm'],
                  "Rh":ratio_dict['Rh']
                  }

        return st, TEROR

    except Exception as e:
        print('\t-- error TERORm creation')
        traceback.print_exc()
        return st, None

def data_gather_and_process(code_set,pre_filt=None,night=False):
    """grab relevant data files for instrument code, process using internal
    functions, return arrays containing filtered waveforms and ratio values
    :type code_set: str
    :param code_set: instrument code set in main
    :type pre_filt: list of two floats
    :param pre_filt: [lower bound, upper bound] filter bounds for wavform plot
    :type night: bool
    :param night: if only night time considered, all files have '_night' suffix
    :rtype st_proc: obspy stream
    :return st_proc: processed stream object for plotting
    :rtype TEORRm: list of numpy arrays
    :return TEORRm: arrays containing filtered waveforms and amp. ratios
    """
    # datapaths
    net,sta,loc,cha,d,year,jday = code_set.split('.')

    fid_path = pathnames()['RDF'] + "{y}/XX/{s}/HH{c}.D".format(y=year,
                                                              s=sta,
                                                              c="{c}")
    inv_path = pathnames()['RDF'] + "DATALESS/XX.RDF.DATALESS"

    # check what combination of files are available, process accordingly
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
                print("\t--too many ({}) traces found...".format(len(st)))
                return None,None
        else:
            print("[stream exists, creating TEROR...]")
            inv = None
            process_check = True
            st = read(path_dict['pickle'])

        # create passband filtered arrays using internal functions
        st_proc, TEORRm = create_modified_TEROR(st,inv,
                                            night=night,
                                            already_preprocessed=process_check)
        if not st_proc:
            print("\t--stream not found...")
            return None,None

        check_bool = check_save(code_set,
                                st=st_proc,
                                TEORRm=TEORRm,
                                night=night)
        if not check_bool:
            return None,None

    else:
        print("[files exist, reading...]")
        st_proc = read(path_dict['pickle'])
        TEORRm = np.load(path_dict['npz'])

    if pre_filt:
        st_proc.filter('bandpass',freqmin=pre_filt[0],freqmax=pre_filt[1])

    return st_proc, TEORRm


# ============================= MAIN PROCESSING ================================
def stacked_process(jday,year='2017'):
    """creating stacked plots by running processing for multiple stations and
    feeding the outputs into plotting script
    data manipualation includes:
    -creating min/max arrays to ease plotting requirements
    -normalizing traces -1 to 1 to allow plotting on the same figure
    -removing values greater than 0.5 (most likely earthquake signals)
    :type jday: int
    :param jday: julian day to process
    :type year: int
    :param year: year to consider, defaults 2017
    """
    # ///////////////////// parameter set \\\\\\\\\\\\\\\\\\\\\\\
    station_list = list(range(1,19))
    night = True
    stop_if_tremor_num_below = 3
    stop_if_stations_above = 16#len(station_list) // 2
    # \\\\\\\\\\\\\\\\\\\\\ parameter set ///////////////////////

    code_set_template = "XX.RD{s:0>2}.10.HH{c}.D.{y}.{d}"
    num_low_tremor_events=0

    # arrays to be filled and passed to plotting function
    y_N_list,y_E_list = [],[]
    TEROR_list,sig_list,tremor_list = [],[],[]
    ano_NE_list = [[],[]]

    # ++ PROCESS EACH STATION AND APPEND TO LISTS
    for station in station_list:
        code_set = code_set_template.format(s=station,c="{c}",d=jday,y=year)
        print('\n\t\t++',code_set)
        # ++ CREATE STREAM OBJECTS AND TEROR ARRAYS
        try:
            st_,TEORRm = data_gather_and_process(code_set,pre_filt=[2,5],
                                                                    night=night)
             # hacky way to avoid overwriting the last entry w/ error
            if not st_:
                print("\t--error no stream object found")
                continue
            st = st_
            st_ = []
        except FileNotFoundError:
            print('\t--file not found')
            continue
        except Exception as e:
            print("\t--error data_gather_and_process()")
            print("="*79)
            traceback.print_exc()
            print("="*79)
            continue

        # ++ COUNT TREMORS AND CREATE SIG ARRAYS
        try:
            tremor_count_dict_ = tremor_counter(TEORRm,
                                                avg_choice="mean",
                                                night=night,
                                                stop_if_tremor_num_below=
                                                    stop_if_tremor_num_below)
            # any false returns means tremor counter exited, tally
            if not tremor_count_dict_:
                num_low_tremor_events += 1
                if stop_if_stations_above and (
                                num_low_tremor_events > stop_if_stations_above):
                    print("Exiting due to low number of station detections")
                    return False
                else:
                    continue
            tremor_count_dict = tremor_count_dict_

        except Exception as e:
            print("\t--error count_tremors_over_sigma()")
            print("="*79)
            traceback.print_exc()
            print("="*79)
            continue

        # ++ SET UP WAVEFORM PLOTTING ARRAYS
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

        sig_list.append(tremor_count_dict)
        TEROR_list.append(TEORRm)


    # ++ CHECK-STOP CONDITIONS
    if not y_N_list or (len(y_N_list) < 3):
        print('\t--not enough stations processed')
        return False

    # ++ OPERATE ON CREATED LISTS
    # RMS of tremor signal for envelope plots, a[a>0.5] == nan
    for N,E in zip(y_N_list,y_E_list):
        horizontal_rms = 0.5 * (N**2 + E**2)
        horizontal_rms /= horizontal_rms.max()
        for data in N,E,horizontal_rms:
            data[abs(data)>0.5] = np.nan
        tremor_list.append(horizontal_rms)

    # ++ TIME - !!! check that all times should be equal eh?
    startNZ,endNZ = convert_UTC_to_local(st)
    t = np.linspace(startNZ.hour,startNZ.hour+12,len(x))

    gridspec_plot(code=code_set,x=t,N=y_N_list,E=y_E_list,
                 TEORRm_list=TEROR_list,
                 sig_list=sig_list,
                 ano_list=ano_NE_list,
                 tremor_list=tremor_list,
                 night=night,
                 show=True,
                 save=False)

    return True


if __name__ == "__main__":
    # already_processed()
    stacked_process(251)
    # for jday in range(215,365):
    #     print('\n==============={}==============='.format(jday))
    #     stacked_process(jday,year='2017')
    # for jday in range(1,150):
    #     print('\n==============={}==============='.format(jday))
    #     stacked_process(jday,year='2018')
