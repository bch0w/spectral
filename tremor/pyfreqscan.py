"""Python port of the modified Frequency Scanning Method (mFSM),
from Katakami et al. 2017 JGR
Used to detect tremors in waveform data by taking amplitude ratio of different
frequency filter bands, removes earthquake data using a simple cross-correlation
with an exponential function. Counts tremor detection using standard deviation
thresholds.
"""
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from obspy import read, read_inventory, Stream
from obspy.signal.cross_correlation import correlate

import sys
sys.path.append("../modules")
from getdata import pathnames

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 0.25
mpl.rcParams['lines.markersize'] = 1.75

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# ============================ PLOTTING FUNCTIONS ==============================
def plot_arrays(st,code_set,TEORRm,sig,show=True):
    """plot 6 streams in 3 subplot figure
    """
    f,(ax1,ax2,ax3) = plt.subplots(3,sharex=True,sharey=False,
                                                    figsize=(9,5),dpi=200)
                                                    
    # divy out arrays
    T,E,O,R,Rm = TEORRm
    sig2,sig3 = sig
                                                    
    # create time axis
    stats = st[0].stats
    start = stats.starttime
    end = stats.endtime
    tracelength = (end-start)/3600
    t = np.linspace(0,tracelength,len(R))
    tRm = np.linspace(0,tracelength,len(Rm))

    # full waveforms
    A1 = ax1.plot(t,st[0].data,color='orange',label=st[0].get_id())
    A2 = ax1.plot(t,st[1].data,color='b',label=st[1].get_id())
    ax1.legend(prop={"size":5})

    # filtered waveforms
    B1 = ax2.plot(t,T,color='r',label="[T]remor Band 2-8Hz")
    B2 = ax2.plot(t,E,color='k',label="[E]arthquake Band 10-20Hz")
    B3 = ax2.plot(t,O,color='g',label="[O]cean Band .5-1Hz")
    ax2.legend(prop={"size":5})

    # amplitude ratio and median value
    # C1 = ax3.plot(t,R,color='k',label='Amplitude [R]atio (T^2/(E*O))')
    C2 = ax3.plot(tRm,Rm,color='k',label='Amplitude [R]atio [m]edian')
    C3 = ax3.scatter(tRm,sig2,c='r',marker='x',zorder=5,label='2-sigma')
    C4 = ax3.scatter(tRm,sig3,c='orange',marker='o',zorder=6,label='3-sigma')
    
    # set axis properties
    for ax in [ax1,ax2,ax3]:
        pretty_grids(ax)
        ax.legend(prop={"size":5})
        
    ax1.set_xlim([t.min(),t.max()])
    st_std = 5 * np.std(st[0].data)
    ax1.set_ylim([-st_std,st_std])
    ax3.set_ylim([Rm.min(),Rm.max()])
    ax2.set_yscale("log")
        
    ax1.set_title("TEROR {}".format(code_set))
    ax1.set_ylabel("velocity (m/s)")
    ax2.set_ylabel("log (velocity (m/s))")
    ax3.set_ylabel("dimensionless")
    ax3.set_xlabel("time (hours from UTC midnight)")

    if show:
        plt.show()

    return f    
    
    
def pretty_grids(input_ax):
    """grid formatting
    """
    input_ax.set_axisbelow(True)
    input_ax.tick_params(which='both',
                         direction='in',
                         top=True,
                         right=True)
    input_ax.minorticks_on()
    input_ax.grid(which='minor',
                    linestyle=':',
                    linewidth='0.5',
                    color='k',
                    alpha=0.25)
    input_ax.grid(which='major',
                    linestyle='-',
                    linewidth='0.5',
                    color='k',
                    alpha=0.15)
    input_ax.ticklabel_format(style='sci',
                            axis='y',
                            scilimits=(0,0))
    
# ========================== SUPPORTING FUNCTIONS ==============================
def check_save(code_set,st=None,TEORRm=None):
    """either check if processing has been run before, or save files to
    specific path. [Check] initiated by default, [save] initiated if the 
    function is given an argument for TEORRm
    :type code_set: str
    :param code_set: instrument code, set in main
    :type st: obspy stream
    :param st: stream containing preprocessed data
    :type TEORRm: list of numpy arrays
    :param TEORRm: list of arrays containing filtered waveform data
    """
    # set up pathing
    net,sta,loc,cha,year,jday = code_set.split('.')
    outpath = pathnames()['data'] + 'TEROR'
    outfile = "{n}.{s}.{l}.{y}.{j}.TEORRm.{f}".format(n=net,
                                                      s=sta,
                                                      l=loc,
                                                      y=year,
                                                      j=jday,
                                                      f='{f}')
    output = os.path.join(outpath,outfile)
    pickle_path = output.format(f='pickle')
    npz_path = output.format(f='npz')
    
    # save arrays
    if TEORRm:
        T,E,O,R,Rm = TEORRm
        st.write(pickle_path,format="PICKLE")
        np.savez(npz_path,T=T,E=E,O=O,R=R,Rm=Rm)
        return True
        
    # check arrays
    else:
        if (os.path.exists(npz_path) and os.path.exists(pickle_path)):
            return {"npz":npz_path,"pickle":pickle_path}
        else:
            return False
        
        
# ========================== PROCESSING FUNCTIONS ==============================
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
    plt.show()
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
            
def data_gather_and_process(code_set):
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
    inv_path = pathnames()['RDF'] + "DATALESS.RDF.XX"
    
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
    :type nighttime: boolean
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
    
    
    

if __name__ == "__main__":
    # ///////////////////// parameter set \\\\\\\\\\\\\\\\\\\\\\\
    code_set_template = "XX.RD06.10.HH{c}.2017.{d}"
    # \\\\\\\\\\\\\\\\\\\\\ parameter set ///////////////////////
    for i in range(213,220):
        code_set = code_set_template.format(c="{c}",d=i)
        print(code_set)
        st,TEORRm,sig = data_gather_and_process(code_set)
        plot_arrays(st,code_set,TEORRm,sig,show=True)
    
        
        
        
        
