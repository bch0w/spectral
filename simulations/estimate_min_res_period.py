"""quick plot script for estimating Minimum Resolveable Period for a coarse
(NGLL5) vs fine (NGLL7) mesh. based on figure 4.7 from Vipul's PhD thesis (2018)
difference between synthetics comes from the GLL points which was set in
constants.h, otherwise both simulations ran on the same parameters
"""
import os
import sys
import glob
sys.path.append('../../modules')
from obspy import read, Stream
from obspy.signal.cross_correlation import correlate, xcorr_max
import numpy as np

from getdata import pathnames
from plotmod import pretty_grids
from synmod import stf_convolve, tshift_halfdur

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1

def get_station_names():
    """get stations in a folder
    """
    # basepath = pathnames()['syns'] + '2018p130600_ngllTest/NGLL5'
    basepath = "/Users/chowbr/Documents/subduction/mseeds/SYN/testing/coarse_2018p130600_ngll/NGLL5"
    checkpath = "/Users/chowbr/Documents/subduction/mseeds/SYN/testing/coarse_2018p130600_ngll/NGLL7"
    filename = "NZ.*.semv.mseed"
    filenames = glob.glob(os.path.join(basepath,filename))

    stalist = []
    for fid in filenames:
        # existpath = os.path.exists(os.path.join(checkpath, os.path.basename(fid)))
        # if existpath:
        stalist.append(fid.split('.')[1])
    stalist = list(set(stalist))


    return stalist

def root_mean_square(u,v):
    """return a single float value of the root mean square for
    """
    numerator = np.linalg.norm(u-v)
    phi = np.sqrt(numerator/np.linalg.norm(u))

    return phi


def estimate_MRP(sta):
    """the function
    """
    # basepath = pathnames()['syns'] + '2018p130600_ngllTest'
    basepath = "/Users/chowbr/Documents/subduction/mseeds/SYN/testing/coarse_2018p130600_ngll"
    objective_filename = "NZ.{s}.{m}XZ.semv.mseed"

    # read in stream objects and mark
    st = Stream()
    for NGLL,marker in zip(['NGLL5','NGLL7'],['B','H']):
        filename = os.path.join(basepath,NGLL,
                                    objective_filename.format(s=sta,m=marker))
        st_tmp = read(filename)
        st_tmp[0].stats.location = NGLL
        st += st_tmp

    # preprocess
    # determine common endtime
    st[1].stats.delta=0.00775
    commonend = min([st[0].stats.endtime,st[1].stats.endtime])
    commonSR = min([st[0].stats.sampling_rate,st[1].stats.sampling_rate])
    for tr in st:
        tr.trim(tr.stats.starttime,commonend)
        tr.resample(commonSR)    
    commonlength = min([len(st[0].data),len(st[1].data)])
    st[0].data = st[0].data[0:commonlength]
    st[1].data = st[1].data[0:commonlength]
    
    # convolve with source time function
    # time_shift,half_duration = tshift_halfdur('2018p130600')
    stn = stf_convolve(st,half_duration=0.7)#,time_shift=time_shift)

    # setup plot
    f,ax = plt.subplots()
    pretty_grids(ax)

    # setup time axis
    t_ngll5 = np.linspace(
                0,stn[0].stats.endtime-stn[0].stats.starttime,len(stn[0].data))
    t_ngll7 = np.linspace(
                0,stn[1].stats.endtime-stn[1].stats.starttime,len(stn[1].data))

    # plot steps of lowpass filters
    step = 0
    for lowpass in np.arange(3,15,0.5):
        st_filter = stn.copy()
        if lowpass != 0:
            st_filter.filter('lowpass',freq=1/lowpass,corners=4)
        for tr in st_filter:
            tr.data /= tr.data.max()
            tr.data += step
        ngll5 = st_filter.select(location='NGLL5')[0].data
        ngll7 = st_filter.select(location='NGLL7')[0].data

        plt.plot(t_ngll5,ngll5,'k',label='NGLL5')
        plt.plot(t_ngll7,ngll7,'r',label='NGLL7')
        # rms = root_mean_square(ngll5,ngll7)
        CC = correlate(ngll5,ngll7,shift=0)

        if lowpass == 0:
            plt.annotate('unfiltered\nCC {:.2f}'.format(CC[0]),xy=(0,step))
        else:
            plt.annotate('lowpass @ {}s\nCC {:.2f}'.format(lowpass,CC[0]),xy=(0,step))

        step += 1.5

    # final plot adjustments
    plt.title('{}'.format(sta)+
                          ' 4KM COARSE 2018p130600 Minimum resolveable period estimation\n'
                                    'NGLL5 (black) vs. NGLL7 (red) Z component')
    plt.gca().axes.get_yaxis().set_ticks([])
    plt.xlabel('Time (s)')
    plt.show()

if __name__ == "__main__":
    stalist = get_station_names()
    for sta in stalist:
        estimate_MRP(sta=sta)
