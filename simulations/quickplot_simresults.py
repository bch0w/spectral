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
from pyatoa.utils.operations.conversions import ascii_to_mseed


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


def process(dir_a, dir_b, min_period=10., max_period=30.,
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
    net, _, comp, ext = os.path.basename(semd).split('.')

    # Remove duplicate stations for search
    stations_fids = np.array([fids, stations])
    unique_stations = np.unique(stations_fids[1, :])

    # Loop through all of the specfem files and read them in
    throwaway_time = UTCDateTime('2000-01-01T00:00:00')
    for sta in unique_stations:
        print(sta)
        tag_template = f"{net}.{sta}.???.{ext}"
        stations_a = glob.glob(os.path.join(dir_a, tag_template))
        stations_b = glob.glob(os.path.join(dir_b, tag_template))

        # Make sure that there are files found
        if not len(stations_a) == len(stations_b):
            continue

        # Make sure the loop finds the correct components
        stations_a.sort()
        stations_b.sort()

        # Convert the data to stream objects
        st_a = Stream()
        st_b = Stream()
        for i in range(len(stations_a)):
            st_a += ascii_to_mseed(stations_a[i], throwaway_time)
            st_b += ascii_to_mseed(stations_b[i], throwaway_time)

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
            ax.plot(t_a, st_a[i].data, color='r', label=st_a[i].get_id())
            ax.plot(t_b, st_b[i].data, color='k', label=st_b[i].get_id())
            ax.legend()
            # Set the title up top
            if i == 0:
                plt.title("{a} (red)\n{b} (black)\n {t0}-{t1}s filter".format(
                    a=os.path.basename(dir_a), b=os.path.basename(dir_b),
                    t0=min_period, t1=max_period)
                )
            # Put the ylabel on the middle plot
            if i == len(st_a) // 2:
                plt.ylabel("amplitude")
            # Put the xlabel on the bottom plot
            if i == len(st_a) - 1:
                plt.xlabel("time [s]")
            # Remove tick labels from all but the last plot
            if i != len(st_a) - 1:
                plt.setp(ax.get_xticklabels(), visible=False)

        if save:
            fid_out = "./{}.png".format(st_a[i].get_id())
            plt.savefig(fid_out, figsize=(11.69, 8.27), dpi=100)
        if show:
            plt.show()


if __name__ == "__main__":
    try:
        base = './'
        dir_a = os.path.join(base, "OUTPUT_FILES_nz_north_coarse_sem_intmesh")
        dir_b = os.path.join(base, "OUTPUT_FILES_srtm30p_139_162_4km_nz_north")
        process(dir_a, dir_b, 10, 30, show=False, save=True)
    except KeyboardInterrupt():
        sys.exit()
