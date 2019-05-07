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


def event_information(pad_s=60*60):
    """
    given an earthquake origin time, produce relevant information
    :return:
    """
    origin_time = UTCDateTime("2017-09-08T04:49:46.0")  # chiapas
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


def beacon_waveforms(number):
    """
    get beacon station waveforms based on station number
    :param number:
    :return:
    """
    start_time, end_time = event_information()
    station_code = "XX.RD{num}.10.HH?.D.{year}.{jday}".format(
        num=number, year=start_time.year, jday=start_time.julday
    )
    net, sta, loc, cha, d, year, jday = station_code.split(".")
    path = "/seis/prj/fwi/bchow/mseeds/BEACON/{year}/XX/{sta}/HH?.D/"
    path = path.format(year=start_time.year, sta=sta)

    st = Stream()
    for fid in glob.glob(os.path.join(path, station_code)):
        st += read(fid)

    st.trim(start_time, end_time)
    inv = read_inventory(
        "/seis/prj/fwi/bchow/mseeds/BEACON/DATALESS/XX.RDF.DATALESS")

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


def plot_components(axes, st):
    """
    plot each stream component on an axis, assuming NEZ
    :param axes:
    :param st:
    :return:
    """
    # assuming all time axes are the same in stream
    for i, component in enumerate(["N", "E", "Z"]):
        tr = st.select(component=component)[0]
        time_axis = np.linspace(0, tr.stats.endtime - tr.stats.starttime,
                                tr.stats.npts
                                )
        axes[i].plot(time_axis, tr.data, label=tr.get_id())
        plt.legend()


def process(t0=1, t1=10):
    """
    main processing
    :return:
    """
    for i in range(1,23,1):
        st_beacon = beacon_waveforms(i)
        st_beacon = preprocess(st_beacon, t0, t1)
        stream_information(st_beacon)
        if st_beacon:
            f, axes = plt.subplots(3, sharex=True)
            plot_components(axes, st_beacon)

            for station in ["BKZ", "PXZ"]:
                st_geonet = geonet_waveforms("NZ.{sta}.10.HH?".format(
                    sta=station)
                )
                st_geonet = preprocess(st_geonet, t0, t1)
                stream_information(st_geonet)
                plot_components(axes, st_geonet)
            plt.show()


if __name__ == "__main__":
    process(t0=10)








