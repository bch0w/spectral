"""
o
Quick script to look at BEACON data for large teleseismic events to check
amplitudes of isntruments in the network
"""

import os
import glob
import pyatoa
import numpy as np
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy.signal.cross_correlation import xcorr_3c
from obspy import UTCDateTime, read, Stream, read_inventory
from pyatoa.utils.tools.process import trimstreams

import matplotlib as mpl
mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['axes.linewidth'] = 2


def change_response(inv):
    from obspy.clients.nrl import NRL
    nrl = NRL()
    resp = nrl.get_response(
        sensor_keys=['Guralp', 'CMG-40T', '60s - 50Hz', '800'],
        datalogger_keys=['Nanometrics', 'Taurus', '40 Vpp (0.4)',
                         'Low (default)', '1 mHz', '100'])
    for net in inv:
        for sta in net:
            for cha in sta:
                cha.response = resp

    return inv


def event_information(start, start_pad, end_pad):
    """
    given an earthquake origin time, produce relevant information
    :return:
    """
    start_time = UTCDateTime(start) - start_pad
    end_time = origin_time + end_pad + start_pad

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


def beacon_waveforms(station, start, end, pad, **kwargs):
    """
    get beacon station waveforms based on station number
    :param number:
    :return:
    """
    path = kwargs.get("path", None)
    inv_path = kwargs.get("inv_path", None)

    code = f"XX.RD{station:0>2}.10.HH?.D.{start.year}.{start.julday:0>3}"
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


def preprocess(st, t0, t1):
    """
    preprocess function with bandpass filter
    :param st:
    :return:
    """
    print(st[0].stats.sampling_rate, end="; ")
    st.decimate(factor=4)
    print(f"decimate to {st[0].stats.sampling_rate}")
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


def plot_components(axes, st, anno='', time_shift=0, color='k', zorder=10):
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

    # assuming all time axes are the same in stream
    for i, component in enumerate(["N", "E", "Z"]):
        tr = st.select(component=component)[0]
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
        axes[i].set_ylabel(f"{component} Amplitude (m/s)")

    axes[1].legend(prop={'size': 6})
    axes[2].set_xlabel("Time since origin time (sec)")
    plt.sca(axes[0])


def test_response_files():
    """
    check amplitude differences in response files
    :return:
    """
    from obspy.clients.nrl import NRL
    nrl = NRL()
    resp_a = nrl.get_response(
        sensor_keys=['Guralp', 'CMG-40T', '60s - 50Hz', '800'],
        datalogger_keys=['Nanometrics', 'Taurus', '40 Vpp (0.4)',
                         'Low (default)', '1 mHz', '100'])
    resp_b = nrl.get_response(
        sensor_keys=['Guralp', 'CMG-40T', '60s - 50Hz', '800'],
        datalogger_keys=['Nanometrics', 'Taurus', '40 Vpp (0.4)',
                         'Low (default)', 'Off', '100'])

    config = pyatoa.Config(
        event_id="2018p130600",
        cfgpaths={
            "waveforms": '/Users/chowbr/Documents/subduction/seismic',
            "responses": '/Users/chowbr/Documents/subduction/seed/RESPONSE'},
        startpad=200, endpad=300
    )
    mgmt = pyatoa.Manager(config)
    rd = mgmt.gatherer.gather_observed("XX.RD10.*.HH?")

    # change response for inv a
    inv_a = mgmt.gatherer.gather_station("XX.RD10.*.HH?")
    for net in inv_a:
        for sta in net:
            for cha in sta:
                cha.response = resp_a
    rd_a = rd.copy()
    rd_a = preprocess(rd_a, inv_a)

    # change response for inv b
    inv_b = mgmt.gatherer.gather_station("XX.RD10.*.HH?")
    for net in inv_b:
        for sta in net:
            for cha in sta:
                cha.response = resp_b
    rd_b = rd.copy()
    rd_b = preprocess(rd_b, inv_b)

    # plot
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    for ax, comp in zip([ax1, ax2, ax3], ["N", "E", "Z"]):
        rd_a_comp = rd_a.select(component=comp)
        rd_b_comp = rd_b.select(component=comp)

        ax.plot(rd_a_comp[0].data, 'r', label="Response A")
        ax.plot(rd_b_comp[0].data, 'k', label="Response B")
    plt.legend()
    plt.show()


def process(starttime, endtime, t0=1, t1=10, normalize=True,
            shift_s=0, pad=50, save=False, show=True, **kwargs):
    """
    main processing
    :return:
    """
    # get the GeoNet data first
    geonet = []
    max_values = {"E": [], "N": [], "Z": []}
    for station in ["BKZ", "PXZ", "TSZ"]:
        st_a = geonet_waveforms(station_code=f"NZ.{station}.10.HH?",
                                start=starttime, end=endtime, pad=pad)
        st_a = preprocess(st_a, t0, t1)

        # Get the max values of each geonet station
        for comp in max_values.keys():
            tr = st_a.select(component=comp)[0]
            max_values[comp].append(tr.data.max())

        geonet.append(st_a)

    # get the Beacon data for each station
    for i in range(1, 23, 1):
        st_b = beacon_waveforms(station=i, start=starttime, end=endtime,
                                pad=pad, **kwargs)
        # If there is no data, continue
        if not st_b:
            print(f"\tNo data for RD{i:0>2}")
            continue
        else:
            print(f"\tRD{i:0>2}")

        st_b = preprocess(st_b, t0, t1)

        # normalize the beacon waveforms by the geonet waveforms
        anno_extra = ''
        if normalize:
            anno_extra = 'normalized'
            for comp in max_values.keys():
                data_ = st_b.select(component=comp)[0].data
                data_max = data_.max()
                data_ /= data_max
                data_ *= max(max_values[comp])
                st_b.select(component=comp)[0].data = data_

        # plot the beacon station
        f, axes = plt.subplots(3, sharex=True)
        plot_components(axes, st_b, color="k", zorder=30, 
                        anno_extra=anno_extra)

        # cross correlate and shift against the main stream so phases match
        for st_, color in zip(geonet, ["r", "b", "orange"]):
            if shift_s:
                st_main = geonet[0].copy()
                st_main.trim(starttime, endtime)
                st_.trim(starttime, endtime)

                shift, value = xcorr_3c(st_main, st_, shift_len=int(
                    st_[0].stats.sampling_rate * shift_s)
                              )
                print(f"{st_[0].get_id()} shift: {shift}s value: "
                      f"{st_[0].stats.sampling_rate/value:.2f}")
                time_shift = shift/st_[0].stats.sampling_rate
            else:
                time_shift = 0
            plot_components(axes, st_, time_shift=time_shift, color=color,
                            zorder=20)

        plt.title(f"{starttime} T=[{t0}, {t1}]s\n"
                  f"{st_b.select(component=comp)[0].get_id()} "
                  f"ratio={data_max/max(max_values[comp]):.2E}")
        if save:
            plt.savefig(f"{starttime.year}_{starttime.month}_rd{i:0>2}.png")
        if show:
            plt.show()

        plt.close()


if __name__ == "__main__":
    # GNS
    path = "/seis/prj/fwi/bchow/data/mseeds/BEACON/{year}/XX/{sta}/HH?.D/"
    inv_path = "/seis/prj/fwi/bchow/data/mseeds/BEACON/DATALESS/beacon.xml"
    # VUW
    path = "/Users/chowbr/Documents/subduction/seismic/mseeds/" \
           "BEACON/{year}/XX/{sta}/HH?.D/"
    inv_path = "/Users/chowbr/Documents/subduction/seismic/mseeds/BEACON/"\
               "DATALESS/beacon.xml"

    origin_times = ["2018-02-18T07:43:48.0",  # 2018p130600 M5.15
                    "2017-09-08T04:49:46.0",  # chiapas, mw8.2
                    "2017-09-19T18:14:48.2",  # central mexico, mw7.1
                    "2018-01-23T09:32:00.0",  # alaska, mw7.9
                    "2018-02-25T17:45:08.6"   # png, mw7.5
                    ]
    for i in range(0,1):
        starttime, endtime = event_information(i, pad=50, end=60*3)
        print(starttime)
        process(starttime, endtime, path=path, normalize=True, 
                inv_path=inv_path, t0=10, t1=30, save=True, show=True)









