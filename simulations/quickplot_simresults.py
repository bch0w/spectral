"""plotting specfem outputs to check the differences
"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from obspy import UTCDateTime
from pyatoa.utils.operations.conversions import ascii_to_mseed


def linespecs():
    import matplotlib as mpl
    mpl.rcParams['font.size'] = 12
    mpl.rcParams['lines.linewidth'] = 1.5
    mpl.rcParams['lines.markersize'] = 1.75
    mpl.rcParams['axes.linewidth'] = 2.0


def pretty_grids(input_ax, scitick=False):
    """make dem grids pretty
    """
    import matplotlib.ticker as ptick
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
                    linestyle=':',
                    linewidth='0.5',
                    color='k',
                    alpha=0.25)
    if scitick:
        input_ax.ticklabel_format(style='sci',
                                axis='y',
                                scilimits=(0,0))


def process():
    """
    if signal processing is required
    """
    base = './'
    folder_1 = "87a78d_old_configure_coarse_mesh"
    folder_2 = "easybuild_19.04_coarse_new_dt"

    output_1 = os.path.join(base,folder_1,'*.sem?')
    semd1 = glob.glob(output_1)
    semd1.sort()

    semd2 = []
    for fid1 in semd1:
        check1 = os.path.basename(fid1)
        check1 = check1.replace('HX', 'BX')
        output_2 = os.path.join(base, folder_2, check1)
        if os.path.exists(output_2):
            semd2.append(output_2)
        else:
            semd1.remove(fid1)
    linespecs()
    for ms1, ms2 in zip(semd1, semd2):
        print(os.path.basename(ms1), os.path.basename(ms2))
        throwaway_time = UTCDateTime('2000-01-01T00:00:00')
        st1 = ascii_to_mseed(ms1, throwaway_time)
        st2 = ascii_to_mseed(ms2, throwaway_time)
        
        t1 = np.linspace(st1[0].stats.time_offset,
                         st1[0].stats.endtime - st1[0].stats.starttime,
                         st1[0].stats.npts)
        t2 = np.linspace(st2[0].stats.time_offset,
                         st2[0].stats.endtime - st2[0].stats.starttime,
                         st2[0].stats.npts) 

        f,ax = plt.subplots(1)
        for st in [st1, st2]:
            st.taper(max_percentage=0.05)
            st.filter('bandpass', freqmin=1/30, freqmax=1/10)
        
        plt.plot(t1, st1[0].data, color='r', label=folder_1)
        plt.plot(t2, st2[0].data, color='k', label=folder_2)
        
        plt.title("{} (red) vs.\n {} (black); 10-30s filter".format(folder_1, folder_2))
        plt.xlabel("Time")
        plt.ylabel("Amplitude")
        pretty_grids(ax)
        
        plt.show()
        sys.exit()


if __name__ == "__main__":
    process()
