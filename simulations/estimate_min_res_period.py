"""
Auick plot script for estimating minimum resolveable period for a coarse
(NGLL5) vs fine (NGLL7) mesh. 

Based on figure 4.7 from Vipul's PhD thesis (2018) difference between synthetics 
comes from the GLL points which was set in constants.h, otherwise both 
simulations ran on the same parameters
"""
import os
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from obspy import read, Stream
from obspy.signal.cross_correlation import correlate, xcorr_max

from pyatoa.utils.tools.io import ascii_to_mseed

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1


def station_names(path_to_syns):
    """
    Station names to be used 
    """
    filenames = glob.glob(os.path.join(basepath, "*.sem?"))

    stalist = []
    for fid in filenames:
        stalist.append(fid.split('.')[1])

    stalist = list(set(stalist))

    return stalist


def estimate_mrp(st_a, st_b, comp="Z", t_min=3, t_max=15):
    """
    estimate the minimum resolvable period between two waveforms based on
    sequential filtering of both traces
    """
    # setup plot
    f,ax = plt.subplots()

    # plot steps of lowpass filters
    step = 1.5
    for lowpass in np.arange(t_min, t_max, 0.5):
        st_a_filt = st_a.copy()
        st_b_filt = st_b.copy()

        # filter
        if lowpass != 0:
            st_a_filt.filter('lowpass', freq=1/lowpass, corners=4)
            st_b_filt.filter('lowpass', freq=1/lowpass, corners=4)

        # normalize all traces, and put them on a step so all plots together
        for tr in st_filter:
            tr.data /= tr.data.max()
            tr.data += step

        data_a = st_a_filt.select(component=comp)[0].data
        data_b = st_b_filt.select(component=comp)[0].data

        plt.plot(t_ngll5,ngll5,'k',label='NGLL5')
        plt.plot(t_ngll7,ngll7,'r',label='NGLL7')

        cc = correlate(data_a, data_b, shift=0)
        shift, value = xcorr_max(cc)

        plt.annotate(f'lowpass @ {lowpass}s\nCC {value:.2f}')

        step += 1.5

    plt.title('Minimum Resolvable Period')
    plt.gca().axes.get_yaxis().set_ticks([])
    plt.xlabel('Time (s)')
    plt.show()


if __name__ == "__main__":
    path_a = ""
    path_b = ""
    stalist = get_station_names(path_a)

    for sta in stalist:
        st_a = ascii_to_mseed(os.path.join(path_a, sta))
        st_b = ascii_to_mseed(os.path.join(path_b, sta))

        estimate_mrp(st_a, st_b)




