"""
Quick script to look at BEACON data for large teleseismic events to check
amplitudes of isntruments in the network
"""

import os
import glob
import json
import pyatoa
import numpy as np
import warnings
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy.signal.cross_correlation import xcorr_3c
from obspy import UTCDateTime, read, Stream, read_inventory
from pyatoa.utils.tools.process import trimstreams

import matplotlib as mpl
mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = .75
mpl.rcParams['axes.linewidth'] = 2

# ignore warnings
warnings.filterwarnings("ignore")


def event_information(start, start_pad, end_pad):
    """
    given an earthquake origin time, produce relevant information
    :return:
    """
    start_time = UTCDateTime(start) - start_pad
    end_time = start_time + end_pad + start_pad

    return start_time, end_time


def geonet_waveforms(station_code, start, end):
    """
    get waveforms from geonet for comparison against BEACON
    :return: 
    """
    net, sta, loc, cha = station_code.split('.')
    c = Client("GEONET")
    st = c.get_waveforms(network=net, station=sta, location=loc, channel=cha,
                         starttime=start, endtime=end, attach_response=True)
    return st


def beacon_waveforms(station_name, start, end, **kwargs):
    """
    get beacon station waveforms based on station number
    :param number:
    :return:
    """
    path = kwargs.get("path", None)
    inv_path = kwargs.get("inv_path", None)

    code = f"XX.{station_name}.10.HH?.D.{start.year}.{start.julday:0>3}"
    net, sta, loc, cha, d, year, jday = code.split(".")

    path = path.format(year=start.year, sta=sta)

    st = Stream()
    for fid in glob.glob(os.path.join(path, code)):
        st += read(fid)
    
    st.trim(start, end)
    inv = read_inventory(inv_path)

    # Will only attach the relevant response
    st.attach_response(inv)

    return st


def preprocess(st_in, t0, t1):
    """
    preprocess function with bandpass filter
    :param st:
    :return:
    """
    st = st_in.copy()

    st.decimate(factor=4)
    st.remove_response(output="VEL", water_level=60,
                       pre_filt=[0.001, 0.005, 45, 50]
                       )
    st.detrend("linear")
    st.detrend("demean")
    st.taper(max_percentage=0.05)

    st.filter("bandpass", freqmin=1/t1, freqmax=1/t0)

    st.detrend("linear")
    st.detrend("demean")
    st.taper(max_percentage=0.05)

    return st


def plot_components(axes, st_in, anno='', time_shift=0, color='k', zorder=10,
                    normalize=True):
    """
    plot each stream component on an axis, assuming NEZ
    :param axes:
    :param st:
    :return:
    """
    def peak_pointer(x, y):
        """
        plot the peak point on the trace
        :param ax:
        :param tr:
        :return:
        """
        peak_y = y.max()
        peak_x = np.where(y == y.max())[0][0]

        return peak_x, peak_y

    # avoid in place edits
    st = st_in.copy()

    # assuming all time axes are the same in stream
    for i, component in enumerate(["N", "E", "Z"]):
        tr = st.select(component=component)[0]
        if normalize:
            tr.data /= tr.data.max()
        time_axis = np.linspace(
            0 + time_shift,
            tr.stats.endtime - tr.stats.starttime + time_shift,
            tr.stats.npts
        )
        axes[i].plot(time_axis, tr.data, label=f"{tr.get_id()} {anno}",
                     color=color, zorder=zorder)
        peak_x, peak_y = peak_pointer(time_axis, tr.data)
        axes[i].scatter(x=time_axis[peak_x], y=peak_y, s=50, c=color,
                        marker="o", edgecolor="k", linewidth=1.5)

        axes[i].grid(True)
        axes[i].set_ylabel(f"{component}")
        axes[i].set_xlim([time_axis[0], time_axis[-1]])

    axes[1].legend(prop={'size': 6})
    axes[2].set_xlabel("Time since origin time (sec)")
    plt.sca(axes[0])


def process(starttime, endtime, t0, t1, plot=False, show=True, save=False,
            normalize=True, **kwargs):
    """
    main processing
    :return:
    """
    # get the GeoNet data first
    geonet = []
    comp_dict = ["N", "E", "Z"]
    amplitude_dict = {}

    for station in ["BKZ", "PXZ", "TSZ"]:
        amplitude_dict[station] = {}

        st_a = geonet_waveforms(station_code=f"NZ.{station}.10.HH?",
                                start=starttime, end=endtime)
        st_a = preprocess(st_a, t0, t1)

        # Get the max values of each geonet station
        for comp in comp_dict:
            tr = st_a.select(component=comp)[0]
            time_ = np.where(tr.data == tr.data.max())[0][0]
            time_ /=  tr.stats.sampling_rate
            amplitude_dict[station][comp] = {"amp": tr.data.max(),
                                             "time": time_}

        geonet.append(st_a)

    # get the Beacon data for each station
    for i in range(1, 23, 1):
        station = f"RD{i:0>2}"
        st_b = beacon_waveforms(station_name=station, start=starttime,
                                end=endtime, **kwargs)
        # If there is no data, continue
        if not st_b:
            print(f"\tNo data for {station}")
            continue
        else:
            print(f"\t{station}")

        amplitude_dict[station] = {}
        st_b = preprocess(st_b, t0, t1)

        for comp in comp_dict:
            tr = st_b.select(component=comp)[0]
            time_ = np.where(tr.data == tr.data.max())[0][0]
            time_ /= tr.stats.sampling_rate
            amplitude_dict[station][comp] = {"amp": tr.data.max(),
                                             "time": time_}
        # Plot the data if necessary
        if plot:
            title = (f"{starttime} T=[{t0}, {t1}]s\n"
                     f"{st_b.select(component=comp)[0].get_id()}"
                     )
            plot_data(geonet=geonet, st_b=st_b, title=title, show=show,
                      save=save, normalize=normalize)

    return amplitude_dict


def plot_data(geonet, st_b, title='', normalize=True, show=True, save=False):
    """
    plot the station data

    :param geonet:
    :param st_b:
    :return:
    """
    # plot the beacon station
    f, axes = plt.subplots(3, sharex=True)
    plot_components(axes, st_b, color="k", zorder=30, normalize=normalize)

    # plot geonet stations
    for color, st_a in zip(["r", "b", "g"], geonet):
        plot_components(axes, st_a, color=color, zorder=20, normalize=normalize)

    starttime = st_b[0].stats.starttime
    plt.title(title)

    if save:
        plt.savefig(
            f"{starttime.year}_{starttime.month:0>2}_{st_b[0].get_id()}.png",
            dpi=100
        )
    if show:
        plt.show()

    plt.close()


def analyze_amp_dict(amp_dict, fid):
    """
    parse through the amplitude dictionary and get amplitude ratios

    :param amp_dict:
    :return:
    """
    with open(fid, "w") as f:
        alreadydone = []
        for sta in amp_dict.keys():
            # Only use geonet stations as baseline comparisons
            if "RD" not in sta:
                continue
            # loop through stations to compare, ignore station weve done
            for sta_compare in amp_dict.keys():
                if (sta == sta_compare) or ("RD" in sta_compare):
                    continue

                combined = sorted(sta + sta_compare)
                jumble = ''.join(combined)
                if jumble in alreadydone:
                    continue
                else:
                    alreadydone.append(''.join(combined))

                print(f"{sta}/{sta_compare} - ", end=' ')
                f.write(f"{sta}/{sta_compare}\t")
                for comp in amp_dict[sta].keys():
                    ratio = (amp_dict[sta][comp]['amp'] /
                             amp_dict[sta_compare][comp]['amp']
                             )
                    print(f"{comp}: {ratio:.2E},", end=' ')
                    f.write(f"{comp}: {ratio:.2E}\t")
                print('')
                f.write("\n")
            print('\n')
            f.write("\n")


if __name__ == "__main__":
    # GNS
    path = "/seis/prj/fwi/bchow/data/mseeds/BEACON/{year}/XX/{sta}/HH?.D/"
    inv_path = "/seis/prj/fwi/bchow/data/mseeds/BEACON/DATALESS/beacon.xml"
    if not os.path.exists(inv_path):
        # VUW
        path = "/Users/chowbr/Documents/subduction/seismic/mseeds/" \
               "BEACON/{year}/XX/{sta}/HH?.D/"
        inv_path = "/Users/chowbr/Documents/subduction/seismic/mseeds/BEACON/"\
                   "DATALESS/beacon.xml"

    origin_times = [("2018-02-18T07:43:48.0", "2018p130600"),  # 2018p130600 M5.15
                    ("2017-09-08T04:49:46.0", "chiapas"),  # chiapas, mw8.2
                    ("2018-01-23T09:32:00.0", "alaska"),  # alaska, mw7.9
                    ("2018-02-25T17:45:08.6", "png"),   # png, mw7.5
                    ("2017-09-19T18:14:48.2", "mexico"),  # central mexico, mw7.1
                    ]
    for event in origin_times[1:]:
        start, end = event_information(start=event[0], start_pad=0,
                                       end_pad=60*60)
        print(start)
        if not os.path.exists(f"./{event[1]}.json"):
            amp_dict = process(start, end, path=path, inv_path=inv_path,
                               plot=True, t0=10, t1=30, save=True, show=False)

            with open(f"{event[1]}.json", "w") as f:
                json.dump(amp_dict, f, indent=4)

        else:
            amp_dict = json.load(open(f"{event[1]}.json"))

        analyze_amp_dict(amp_dict)








