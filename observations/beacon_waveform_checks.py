"""
Quick script to look at BEACON data for large teleseismic events to check
amplitudes of isntruments in the network
"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy import UTCDateTime, read, Stream, read_inventory

import matplotlib as mpl
mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['axes.linewidth'] = 2


def event_information(i=0, pad_s=60*60):
    """
    given an earthquake origin time, produce relevant information
    :return:
    """
    origin_times = ["2017-09-08T04:49:46.0",  # chiapas, mw8.2
                    "2017-09-19T18:14:48.2",  # central mexico, mw7.1
                    "2018-01-23T09:32:00.0",  # alaska, mw7.9
                    "2018-02-25T17:45:08.6"   # png, mw7.5
                    ]
    origin_time = UTCDateTime(origin_times[i])
    end_time = origin_time + pad_s

    return origin_time, end_time


def geonet_waveforms(station_code):
    """
    get waveforms from geonet for comparison against BEACON
    :return: 
    """
    net, sta, loc, cha = station_code.split('.')
    c = Client("GEONET")
    start_time, end_time = event_information()
    st = c.get_waveforms(network=net, station=sta, location=loc, channel=cha,
                         starttime=start_time, endtime=end_time,
                         attach_response=True
                         )
    return st


def beacon_waveforms(number, path, inv_path):
    """
    get beacon station waveforms based on station number
    :param number:
    :return:
    """
    start_time, end_time = event_information(i=1)
    station_code = "XX.RD{num}.10.HH?.D.{year}.{jday}".format(
        num=number, year=start_time.year, jday=start_time.julday
    )
    net, sta, loc, cha, d, year, jday = station_code.split(".")
    path = path.format(year=start_time.year, sta=sta)

    st = Stream()
    for fid in glob.glob(os.path.join(path, station_code)):
        st += read(fid)

    st.trim(start_time, end_time)
    inv = read_inventory(inv_path)

    st.attach_response(inv)

    return st


def preprocess(st, t0, t1):
    """
    preprocess function with bandpass filter
    :param st:
    :return:
    """
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


def stream_information(st):
    """
    print out relevant stream information
    :param st:
    :return:
    """
    for tr in st:
        timing = np.where(tr.data == tr.max())[0][0] / tr.stats.sampling_rate
        print("Channel Code: {id}\n"
              "Starttime: {st}\n"
              "Endtime: {en}\n"
              "Max Amplitude: {amp}\n"
              "Timing of Max: {tim} s\n".format(id=tr.get_id(),
                                                st=tr.stats.starttime,
                                                en=tr.stats.endtime,
                                                amp=tr.max(),
                                                tim=timing)
              )


def plot_components(axes, st, c):
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
        time_axis = np.linspace(0, tr.stats.endtime - tr.stats.starttime,
                                tr.stats.npts
                                )
        axes[i].plot(time_axis, tr.data, label=tr.get_id(), color=c)
        peak_x, peak_y = peak_pointer(time_axis, tr.data)
        axes[i].plot(
            time_axis[peak_x], peak_y, '{}o'.format(c), markersize=4)

        axes[i].grid(True)
        axes[i].set_ylabel("{} Amplitude (m/s)".format(component))

    axes[1].legend()
    axes[2].set_xlabel("Time since origin time (sec)")
    plt.sca(axes[0])


def process(path, inv_path, t0=1, t1=10):
    """
    main processing
    :return:
    """
    for i in range(1, 23, 1):
        st_beacon = beacon_waveforms(i, path, inv_path)
        st_beacon = preprocess(st_beacon, t0, t1)
        stream_information(st_beacon)
        if st_beacon:
            f, axes = plt.subplots(3, sharex=True)
            plot_components(axes, st_beacon, c="b")

            for station, color in zip(["BKZ", "PXZ"], ["r", "k"]):
                st_geonet = geonet_waveforms("NZ.{sta}.10.HH?".format(
                    sta=station)
                )
                st_geonet = preprocess(st_geonet, t0, t1)
                stream_information(st_geonet)
                plot_components(axes, st_geonet, c=color)
            plt.title("{} T=[{},{}]s".format(event_information()[0], t0, t1))
            plt.show()


if __name__ == "__main__":
    pathing = "VUW"
    if pathing == "GNS":
        path = "/seis/prj/fwi/bchow/mseeds/BEACON/{year}/XX/{sta}/HH?.D/"
        inv_path = "/seis/prj/fwi/bchow/mseeds/BEACON/DATALESS/XX.RDF.DATALESS"
    elif pathing == "VUW":
        path = "/Users/chowbr/Documents/subduction/mseeds/" \
               "BEACON/{year}/XX/{sta}/HH?.D/"
        inv_path = "/Users/chowbr/Documents/subduction/mseeds/" \
                   "BEACON/DATALESS/XX.RDF.DATALESS"

    process(path, inv_path, t0=10, t1=30)








