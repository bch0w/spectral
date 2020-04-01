"""
Quick plot script for estimating minimum resolveable period for a coarse
(NGLL5) vs fine (NGLL7) mesh. This gives a rough idea about what the minimum 
synthetic period that can be used is, gives a good idea of what bandpass
the mesh will work for

Based on figure 4.7 from Vipul Silwahl's PhD thesis (2018) difference between 
synthetics comes from the GLL points which was set in constants.h, otherwise
simulations should be run with the same parameters

Cross correlation used to determine the similarity between traces
"""
import os
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from obspy import read, Stream
from obspy.signal.cross_correlation import correlate, xcorr_max

from pyatoa.utils.tools.io import read_ascii
from pyatoa.utils.visuals.plot_tools import pretty_grids

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = .6


def station_names(path_to_syns):
    """
    Station names to be used 
    """
    filenames = glob.glob(os.path.join(path_to_syns, "*.sem?"))

    stalist = []
    for fid in filenames:
        stalist.append(os.path.basename(fid).split('.')[1])

    stalist = list(set(stalist))
    stalist.sort()

    return stalist


def estimate_mrp(st_a, st_b, t_min=3, t_max=15, step_count=1.1, resample=False,
                 ccorr=False):
    """
    estimate the minimum resolvable period between two waveforms based on
    sequential filtering of both traces
    """
    # setup plot
    f,ax = plt.subplots(figsize=(5,8), dpi=200)
    pretty_grids(ax) 

    # plot the raw data first
    plt.plot(st_a[0].times(), st_a[0].data/st_a[0].data.max(), 'k', 
             label='ngll5')
    plt.plot(st_b[0].times(), st_b[0].data/st_b[0].data.max(), 'r', 
             label='ngll7')
    plt.annotate('raw data', xy=(0, 0), fontsize=4)

    # plot steps of lowpass filters
    step = step_count
    for i, lowpass in enumerate(np.arange(t_min, t_max, 0.5)):
        st_a_filt = st_a.copy()
        st_b_filt = st_b.copy()

        # filter
        if lowpass != 0:
            st_a_filt.filter('lowpass', freq=1/lowpass, corners=4)
            st_b_filt.filter('lowpass', freq=1/lowpass, corners=4)
        
        # normalize all traces, and put them on a step so all plots together
        for st in [st_a_filt, st_b_filt]: 
            for tr in st:
                tr.data /= tr.data.max()
                tr.data += step
        
        # find common sampling rate, resample higher sampling rate data
        if resample:
            sr_a = st_a_filt[0].stats.sampling_rate
            sr_b = st_b_filt[0].stats.sampling_rate
            if sr_a > sr_b:
                st_a_filt.resample(sr_b)
            elif sr_a < sr_b:
                st_b_filt.resample(sr_a)
       
        # assume that there is only one trace in the stream
        plt.plot(st_a_filt[0].times(), st_a_filt[0].data, 'k')
        plt.plot(st_b_filt[0].times(), st_b_filt[0].data, 'r')

        if ccorr and resample:
            cc = correlate(st_a_filt[0].data, st_b_filt[0].data, shift=0)
            shift, value = xcorr_max(cc)

            plt.annotate(f'lowpass @ {lowpass}s\nCC {value:.2f}', xy=(0, step),
                         fontsize=4)

        step += step_count

    # Plot accessories
    plt.title(f'Minimum Resolvable Period')
    plt.grid()
    plt.legend(loc="upper right")
    plt.xlabel('Time (s)')
    plt.ylabel(f"Filter bands {t_min}-{t_max}s")

    plt.gca().axes.get_yaxis().set_ticks([])
    plt.savefig(f'mrp_{st_a[0].get_id()}.png')
    plt.close()


if __name__ == "__main__":
    path_a = "./ngll5"
    path_b = "./ngll7"
    comp = "Z"
    stalist = station_names(path_a)

    for sta in stalist:
        fid_template = f"??.{sta}.??{comp}.sem?"

        fid_a = glob.glob(os.path.join(path_a, fid_template))
        fid_b = glob.glob(os.path.join(path_b, fid_template))
        if fid_a and fid_b:
            print(sta)
            st_a = read_ascii(fid_a[0])
            st_b = read_ascii(fid_b[0])

            estimate_mrp(st_a, st_b, t_min=1, resample=True, ccorr=True)




