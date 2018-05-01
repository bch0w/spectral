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
    """
    # create UTC aware datetime objects
    startUTC = st[0].stats.starttime.datetime.replace(tzinfo=pytz.utc)
    endUTC = st[0].stats.endtime.datetime.replace(tzinfo=pytz.utc)
    startLOC = startUTC.astimezone(pytz.timezone(local_timezone))
    endLOC = endUTC.astimezone(pytz.timezone(local_timezone))

    return startLOC, endLOC

def average_array(data,number_of_samples):
    """take a numpy array and return a shortened version which contains averaged
    values for every number of samples specified. the end point may contain a
    disproportionate number of samples to make up for the indivisibility of
    the number of samples by the chosen division rate

    np.nanmean throws a runtime warning for all-nan slices, and returns nan
    """
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    new_data,std_list = [],[]
    for i in range(0,len(data)-number_of_samples,number_of_samples):
        j = i + number_of_samples
        if j > len(data)-number_of_samples:
            j = len(data) - 1
        new_data.append(np.nanmean(data[i:j]))

    return np.array(new_data)


def create_horizontal_data(st,bounds):
    """filter north and east component data and combine into horizontal comp
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

def mean_std_creator(night=False):
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
    outpath = os.path.join(filepath,'averages{N}_{C}.npz'.format(
                                                                N=N,
                                                                C=num_of_files))
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
    outpath_wild = os.path.join(filepath,'averages{N}_{C}.npz'.format(
                                                                    c=choice,
                                                                    N=N,
                                                                    C="*"))
    files_to_delete = glob.glob(outpath_wild)
    files_to_delete.remove(outpath)
    if files_to_delete:
        for fid in files_to_delete:
            os.remove(fid)

    return sig_dict

# =========================== PROCESSING FUNCTIONS =============================

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
    st_copy = st_raw.copy()
    T0 = time.time()
    # trim traces to nz local night time, NZ ~= UTC+12
    if night:
        year = st_raw[0].stats.starttime.year
        month = st_raw[0].stats.starttime.month
        day = st_raw[0].stats.starttime.day

        sunset_in_nz = UTCDateTime(year,month,day,6,0,0)
        sunrise_in_nz = UTCDateTime(year,month,day,18,0,0)

        # sunset_in_nz = st_raw[0].stats.starttime + 6*3600
        # sunrise_in_nz = st_raw[0].stats.starttime + 18*3600
        st_raw.trim(sunset_in_nz,sunrise_in_nz)

    # data quality checks
    if not st_raw:
        print("streams fall outside trim")
        return False
    elif st_raw[0].stats.npts != st_raw[1].stats.npts:
        print("streams have inconsistent npts {0} {1}".format(
                                                        st_raw[0].stats.npts,
                                                        st_raw[1].stats.npts))
        return False

    decimate_by = int(st_raw[0].stats.sampling_rate//resample)
    st_raw.decimate(decimate_by)
    st_raw.detrend("demean") # not in original code
    st_raw.detrend("linear") # not in original code
    st_raw.taper(max_percentage=0.05)

    st_raw.attach_response(inv)
    pre_filt = [.001,.01,50,55]
    st_raw.remove_response(output="VEL",
                          pre_filt=pre_filt, # not in original code
                          water_level=water_level, # not in original code
                          plot=False)

    print(round(time.time()-T0,2))
    return st_raw

def detect_earthquakes(tremor_horizontal,sampling_rate,corr_criteria=0.7):
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
    sampling_rate_min = sampling_rate * 60
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
    nan_fill = np.nan*(np.ones(sampling_rate_min))
    nan_fill_ext = np.nan*(np.ones(sampling_rate_one_one_half_min))

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
                quakearray = np.append(quakearray,nan_fill)
            else:
                quakearray_new = quakearray[:S0-sampling_rate_half_min]
                quakearray = np.append(quakearray_new,nan_fill_ext)
        else:
            quakearray = np.append(quakearray,tremor_snippet)

    print("{} earthquakes detected".format(quakecount))

    return quakearray

def tremor_counter(array,choice,avg_choice='mean',night=False):
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
    sig_dict = mean_std_creator(night)
    mean_ = sig_dict["{}_mean".format(choice)]
    median_ = sig_dict["{}_median".format(choice)]
    one_sigma = sig_dict["{}_sigma".format(choice)]

    avg_val = mean_
    if avg_choice == 'median':
        avg_val = median_

    two_sigma = avg_val + one_sigma * 2
    three_sigma = avg_val + one_sigma * 3

    R_2sig = np.copy(array)
    R_3sig = np.copy(array)

    # count tremor activity by fulfilling threshold criteria for 2sig and 3sig
    for i,section in enumerate(array):
        if section < two_sigma:
            R_2sig[i] = 0
            R_3sig[i] = 0
        elif section < three_sigma:
            R_3sig[i] = 0

    print("[tremor_counter {c}: 2-sig:{two} / 3-sig:{three}]".format(c=choice,
                                                    two=len(R_2sig[R_2sig>0]),
                                                    three=len(R_3sig[R_3sig>0])
                                                    ))

    sigma_dict = {"2-sigma":R_2sig,
                  "3-sigma":R_3sig}

    return sigma_dict

def count_tremors_over_sigma(TEORRm,avg_choice="median",night=False,
                                               stop_if_tremor_num_below = 3):
    """count number of tremors in a given array and return values which cross
    a global 2-sigma and 3-sigma threshold, which signifies a tremor 'detection'
    station specific function. Currently only returns 2-sigma values
    :return tremor_count_dict: dictionary filled with arrays where ratio arrays
    surpass 2-sigma over the mean value
    """
    tremor_count_dict = {}
    for choice in ['Rs','Rm','Rh']:
        sgtmp = tremor_counter(TEORRm[choice],choice=choice,
                                        avg_choice=avg_choice,night=night)
        two_sigma = sgtmp["2-sigma"]

        # check-stop if low-tremor detection threshold met on Rs
        if stop_if_tremor_num_below and choice=='Rm':
            if len(two_sigma[two_sigma>0]) < stop_if_tremor_num_below:
                print("\t--tremor number below {n} ".format(
                                                n=stop_if_tremor_num_below))
                return None

        # nan out 0's to avoid overplotting 0 values
        two_sigma[two_sigma==0] = np.nan
        tremor_count_dict[choice] = two_sigma

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
        tremor_horizontal = create_horizontal_data(st,tremor_band_short)
        earthquake_horizontal = create_horizontal_data(st,earthquake_band_short)
        surface_wave_horizontal = create_horizontal_data(st,surface_wave_band)

        # simple earthquake detection and removal
        tremor_horizontal = detect_earthquakes(tremor_horizontal,sampling_rate,
                                                            corr_criteria=0.65)

        # create ratio arrays by certain time windows,
        # honor earthquake detection in all ratios, replace any nan's with 0's
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

        # zero out any values of Rm > (Rs.mean + 3*Rs.1-sigma)
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
                print("{} traces found, skipping".format(len(st)))
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
    """
    # ///////////////////// parameter set \\\\\\\\\\\\\\\\\\\\\\\
    station_list = list(range(0,19))
    nightonly = True
    stop_if_tremor_num_below = 3
    stop_if_stations_above = len(station_list) // 2
    # \\\\\\\\\\\\\\\\\\\\\ parameter set ///////////////////////

    # accumulate, process all data and place in arrays for later plotting
    code_set_template = "XX.RD{s:0>2}.10.HH{c}.D.{y}.{d}"
    num_low_tremor_events=0

    y_N_list,y_E_list,TEROR_list,sig_list,tremor_list = [],[],[],[],[]
    ano_NE_list = [[],[]]

    # ++ PROCESS EACH STATION AND APPEND TO LISTS
    for station in station_list:
        code_set = code_set_template.format(s=station,c="{c}",d=jday,y=year)
        print('\n',code_set)
        # ++ CREATE STREAM OBJECTS AND TEROR ARRAYS
        try:
            st_,TEORRm = data_gather_and_process(code_set,pre_filt=[2,5],
                                                                night=nightonly)
             # hacky way to avoid overwriting the last entry w/ error
            if not st_:
                print("\t--error no stream object found")
                continue
            st = st_
            st_ = []
        except FileNotFoundError:
            print('\t--nonexistent file')
            continue
        except Exception as e:
            print("\t--error data_gather_and_process")
            print("="*79)
            traceback.print_exc()
            print("="*79)
            continue

        # ++ COUNT TREMORS AND CREATE SIG ARRAYS
        try:
            tremor_count_dict_ = count_tremors_over_sigma(TEORRm,
                                            avg_choice="median",night=False,
                                            stop_if_tremor_num_below=
                                                    stop_if_tremor_num_below)
            if not tremor_count_dict_:
                num_low_tremor_events += 1
                if stop_if_stations_above and (
                                num_low_tremor_events > stop_if_stations_above):
                    print("Exiting due to low number of station detections")
                    return False
            if not tremor_count_dict_:
                continue
            tremor_count_dict = tremor_count_dict_
        except Exception as e:
            print("[error count_tremors_over_sigma]")
            print("="*79)
            traceback.print_exc()
            print("="*79)
            continue

        # ++ SET UP WAVEFORM PLOTTING ARRAYS
        for comp in ["N","E"]:
            print(len(st[0]))
            x,y = create_min_max(st.select(component=comp)[0])
            y_E_max = y.max()
            if comp == "N":
                y_N_max = y.max()
            y/=y.max()
            y_N_list.append(y) if comp == "N" else y_E_list.append(y)



        sig_list.append(tremor_count_dict)
        TEROR_list.append(TEORRm)

        # create annotations from max amplitude values
        for i,max_ in enumerate([y_N_max,y_E_max]):
            ano_NE_list[i].append("{s} {a}um/s".format(s=code_set.split('.')[1],
                                                        a=round(max_ * 1E6,2)
                                                        ))
    # check stop conditions
    if not y_N_list or (len(y_N_list) < 3):
        print('\t--not enough stations processed')
        return False

    # ++ OPERATE ON CREATED LISTS
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

    gridspec_plot(code=code_set,x=t,N=y_N_list,E=y_E_list,
                 TEORRm_list=TEROR_list,
                 sig_list=sig_list,
                 ano_list=ano_NE_list,
                 tremor_list=tremor_list,
                 night=nightonly,
                 show=True,
                 save=False)

    return True



if __name__ == "__main__":
    # already_processed()
    stacked_process(251)
    # for jday in range(199,365):
    #     print('\n==============={}==============='.format(jday))
    #     stacked_process(jday,year='2017')
    # for jday in range(1,63):
    #     print('\n==============={}==============='.format(jday))
    #     stacked_process(jday,year='2018')
