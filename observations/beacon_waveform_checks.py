"""
Quick script to look at BEACON data for large teleseismic events to check
amplitudes of isntruments in the network
"""
import os
import glob
import pyatoa
import numpy as np
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy import UTCDateTime, read, Stream, read_inventory

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


def compare_geonet():
    """
    Comparing BEACON stations to GeoNet sites nearby to check waveform fits.
    :return:
    """
    config = pyatoa.Config(event_id="2017p795065",
                           paths_to_waveforms=[
                               '/Users/chowbr/Documents/subduction/seismic'],
                           paths_to_responses=[
                               '/Users/chowbr/Documents/subduction/seed/RESPONSE'],
                           startpad=200, endpad=300
                           )
    mgmt = pyatoa.Manager(config)
    # NZ.BFZ
    bfz = mgmt.gatherer.gather_observed(station_code="NZ.TSZ.*.HH?")
    bfz_inv = mgmt.gatherer.gather_station(
        station_code="NZ.TSZ.*.HH?")
    bfz = preprocess(bfz, bfz_inv)
    # NZ.BKZ
    bkz = mgmt.gatherer.gather_observed(station_code="NZ.WAZ.*.HH?")
    bkz_inv = mgmt.gatherer.gather_station(
        station_code="NZ.WAZ.*.HH?")
    bkz = preprocess(bkz, bkz_inv)
    # NZ.PXZ
    pxz = mgmt.gatherer.gather_observed(station_code="NZ.PXZ.*.HH?")
    pxz_inv = mgmt.gatherer.gather_station(
        station_code="NZ.PXZ.*.HH?")
    pxz = preprocess(pxz, pxz_inv)

    for i in range(1, 20):
        try:
            f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
            rd = mgmt.gatherer.gather_observed(
                "XX.RD{:0>2}.*.HH?".format(i))
            rd_inv = mgmt.gatherer.gather_station(
                "XX.RD{:0>2}.*.HH?".format(i))
            rd = preprocess(rd, rd_inv)

            for ax, comp in zip([ax1, ax2, ax3], ["N", "E", "Z"]):
                bfz_comp = bfz.select(component=comp)
                bkz_comp = bkz.select(component=comp)
                pxz_comp = pxz.select(component=comp)
                rd_comp = rd.select(component=comp)

                ax.plot(bfz_comp[0].data, 'k', label="TSZ")
                ax.plot(bkz_comp[0].data, 'k--', label="WAZ")
                ax.plot(pxz_comp[0].data, 'k:', label="PXZ")
                ax.plot(rd_comp[0].data, 'r', label="RD??")

                ax.set_ylabel(comp)
                ax.set_title("RD{:0>2}".format(i))
                plt.legend()
            print("BFZ {bfz}\nBKZ {bkz}\nPXZ {pxz}\nRD {rd}".format(
                bfz=bfz.max(), bkz=bkz.max(), pxz=pxz.max(),
                rd=rd.max())
            )
            plt.show()
        except Exception as e:
            plt.close()
            continue

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

    config = pyatoa.Config(event_id="2018p130600",
                           paths_to_waveforms=[
                               '/Users/chowbr/Documents/subduction/seismic'],
                           paths_to_responses=[
                               '/Users/chowbr/Documents/subduction/seed/RESPONSE'],
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









