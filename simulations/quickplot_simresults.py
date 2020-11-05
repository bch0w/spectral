"""
Quickly plotting specfem outputs to check the differences in sythetics
"""
import os
import sys
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from obspy import UTCDateTime, Stream
from obspy.signal.cross_correlation import correlate, xcorr_max
from pyatoa import read_sem


def linespecs():
    """
    Set line specifications
    :return:
    """
    mpl.rcParams['font.size'] = 10
    mpl.rcParams['lines.linewidth'] = 1.5
    mpl.rcParams['lines.markersize'] = 1.75
    mpl.rcParams['axes.linewidth'] = 2.0


def pretty_grids(input_ax, scitick=False):
    """
    Control the look of the underlying axes and grid
    """
    input_ax.set_axisbelow(True)
    input_ax.tick_params(which='both',
                         direction='in',
                         top=True,
                         right=True)
    input_ax.minorticks_on()
    input_ax.grid(which='minor', linestyle=':', linewidth='0.5', color='k',
                  alpha=0.25)
    input_ax.grid(which='major',linestyle=':', linewidth='0.5', color='k',
                  alpha=0.25)
    if scitick:
        input_ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))


def parse_by_component(dir_a, dir_b):
    """

    :param dir_a:
    :param dir_b:
    :return:
    """
    output_a = os.path.join(dir_a, '*.sem?')
    semd_a = glob.glob(output_a)
    semd_a.sort()

    # Get all the files from directory B and match them to directory A
    semd_b = []
    for fid_a in semd_a:
        check_a = os.path.basename(fid_a)
        # Incase the sampling rates are different
        check_a = check_a.replace('HX', 'BX')

        output_b = os.path.join(dir_b, check_a)
        if os.path.exists(output_b):
            semd_b.append(output_b)
        else:
            semd_a.remove(fid_a)


def peak_amplitudes(tr, t, ax, c, plot=False):
    """
    Plot the location of the peak amplitudes and their time difference
    :param st_a:
    :param st_b:
    :param ax:
    :return:
    """
    peak_amp = tr.data.max()
    time_peak = t[np.where(tr.data == tr.data.max())[0]]

    if plot:
        ax.axvline(x=time_peak, ymin=0, ymax=1, c=c, linestyle="--")

    return time_peak, peak_amp


def process(dir_a, dir_b, min_period=10., max_period=30., peak_amplitude=False,
            cross_correlate=True, color_a="darkorange", color_b="mediumblue",
            show=True, save=False):
    """
    Read and process the data
    """
    linespecs()

    # Get all the specfem files from directory A, sort them alphabetically
    semd_files_a = glob.glob(os.path.join(dir_a, '*.sem?'))
    fids, stations = [], []
    for semd in semd_files_a:
        fids.append(os.path.basename(semd))
        stations.append(os.path.basename(semd).split('.')[1])

    # Set the tag template for glob to search for stations
    _, _, _, ext = os.path.basename(semd).split('.')

    # Remove duplicate stations for search
    stations_fids = np.array([fids, stations])
    unique_stations = np.unique(stations_fids[1, :])

    # Loop through all of the specfem files and read them in
    throwaway_time = UTCDateTime('2000-01-01T00:00:00')
    for sta in unique_stations:
        print(sta, end=" ")
        tag_template = f"??.{sta}.???.{ext}"
        stations_a = glob.glob(os.path.join(dir_a, tag_template))
        stations_b = glob.glob(os.path.join(dir_b, tag_template))

        # Make sure that there are files found
        if not (len(stations_a) == len(stations_b)):
            print("skipped")
            continue

        # Make sure the loop finds the correct components
        stations_a.sort()
        stations_b.sort()

        # Convert the data to stream objects
        st_a = Stream()
        st_b = Stream()
        for i in range(len(stations_a)):
            st_a += read_sem(stations_a[i], throwaway_time)
            st_b += read_sem(stations_b[i], throwaway_time)
        
        if not st_a or not st_b: 
            print("skipped")
            continue
        
        # Create common time axes
        t_a = np.linspace(st_a[0].stats.time_offset,
                          st_a[0].stats.endtime - st_a[0].stats.starttime,
                          st_a[0].stats.npts)
        t_b = np.linspace(st_b[0].stats.time_offset,
                          st_b[0].stats.endtime - st_b[0].stats.starttime,
                          st_b[0].stats.npts)

        # Filter the data and taper to remove any spurious signals
        for st in [st_a, st_b]:
            st.taper(max_percentage=0.05)
            if min_period:
                st.filter('bandpass', freqmin=1/max_period,
                          freqmax=1/min_period)
                st.taper(max_percentage=0.05)

        # Plot each component on a Grid
        gs = gridspec.GridSpec(3, 1, hspace=0)
        for i in range(len(st_a)):
            ax = plt.subplot(gs[i])
            pretty_grids(ax, scitick=True)
            anno = ""

            # Plot the peak amplitude differences, and time differences
            time_peak_a, peak_a = peak_amplitudes(st_a[i], t_a, ax, c=color_a,
                                                  plot=peak_amplitude)
            time_peak_b, peak_b = peak_amplitudes(st_b[i], t_b, ax, c=color_b,
                                                  plot=peak_amplitude)
            if peak_amplitude:
                anno = "delta_amp={:.2E}\ndelta_t={:.1f}s\n".format(
                    float(peak_a - peak_b), float(time_peak_a - time_peak_b))

            # Cross correlate the two traces and annotate the cc information
            if cross_correlate:
                common_sr = min([st_a[i].stats.sampling_rate, 
                                 st_b[i].stats.sampling_rate])
                tr_a = st_a[i].copy()
                tr_b = st_b[i].copy()
                for tr in [tr_a, tr_b]:
                    tr.resample(common_sr) 
                
                cc = correlate(tr_a.data, tr_b.data, shift=int(common_sr * 50), 
                               domain="freq")

                f_shift, value = xcorr_max(cc)
                t_shift = f_shift / common_sr
                anno += f"cc={value:.2f}\ntshift={t_shift:.2f}s"
              
            ax.annotate(s=anno, 
                        xy=(t_a[int(len(t_a)/2)], -1 * st_a[i].data.max()/2), 
                        fontsize=8, 
                        bbox=dict(fc="w", boxstyle="round", alpha=0.5))
                
            # Plot the two traces
            ax.plot(t_a, st_a[i].data, color=color_a)
            ax.plot(t_b, st_b[i].data, color=color_b)
            # plt.xlim([0, 100])

            # Set the title with important information
            if i == 0:
                plt.title(f"2018p130600 {st_a[i].get_id()}\n"
                          f"{os.path.basename(dir_a)} ({color_a}) / "
                          f"{os.path.basename(dir_b)} ({color_b})\n"
                          f"{min_period} - {max_period}s")
            # Put the ylabel on the middle plot
            if i == len(st_a) // 2:
                plt.ylabel(f"amplitude [m]\n{st_a[i].get_id().split('.')[-1]}")
            else:
                plt.ylabel(f"{st_a[i].get_id().split('.')[-1]}")
            # Put the xlabel on the bottom plot
            if i == len(st_a) - 1:
                plt.xlabel("time [s]")
            # Remove tick labels from all but the last plot
            if i != len(st_a) - 1:
                plt.setp(ax.get_xticklabels(), visible=False)

        if save:
            fid_out = "./figures/{}.png".format(st_a[i].get_id())
            plt.savefig(fid_out, figsize=(11.69, 8.27), dpi=100)
        if show:
            plt.show()
        
        print("")


if __name__ == "__main__":
    base = './'
    dir_a = os.path.join(base, "SEMD_HIRES")
    dir_b = os.path.join(base, "SEMD_SUPERHIRES")
    if not os.path.exists(os.path.join(base, "figures")):
        os.makedirs(os.path.join(base, "figures"))

    process(dir_a, dir_b, 8, 30, peak_amplitude=False, 
            cross_correlate=True, show=False, save=True)
