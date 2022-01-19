"""
Check long period GeoNet waveforms to see if there are any amplitude 
reversals for stations. Longer period waveforms should be relatively in phase
"""
import os
import glob
import pprint
import numpy as np
from scipy import signal
from obspy.clients.fdsn import Client
from obspy import read, read_inventory, UTCDateTime, Stream
from obspy.clients.fdsn.header import FDSNException

import matplotlib.pyplot as plt
import matplotlib as mpl


c = Client("GEONET")

class DataGapError(Exception):
    """
    Error resulting from data gaps
    """
    pass


def geonet_waveforms(station_code, start, end):
    """
    Get waveforms from geonet for comparison against BEACON
    """
    print(f"GeoNet waveform {station_code}", end="... ")
    net, sta, loc, cha = station_code.split('.')

    path_to = "./data"
    fid = os.path.join(path_to,
                       f"{station_code}.D.{start.year}.{start.julday:0>3}")
    inv_fid = os.path.join(path_to, f"RESP.{net}.{sta}.{loc}.{cha}")

    # Get Data
    if os.path.exists(fid) and os.path.exists(inv_fid):
        print("internal")
        st = read(fid)
        inv = read_inventory(inv_fid)
    else:
        print("external", end="... ")

        try:
            st = c.get_waveforms(network=net, station=sta, location=loc,
                                 channel=cha, starttime=start, endtime=end,
                                 attach_response=False)
        except FDSNException:
            raise FileNotFoundError(f"GeoNet waveforms could not be fetched")
        try:
            inv = c.get_stations(network=net, station=sta, location=loc,
                                 channel=cha, starttime=start,
                                 endtime=end, level="response")
        except FDSNException:
            raise FileNotFoundError(f"GeoNet inventory could not be fetched")

        print("writing")
        st.write(fid, format="mseed")
        inv.write(inv_fid, format="stationxml")

    if not st:
        raise FileNotFoundError(f"GeoNet data not found for {station_code}")
    if len(st) > 1:
        raise Exception("Non continuous data")

    st.trim(start, end)

    return st, inv


def preprocess(st_in, inv=None, sampling_rate=None, response=True):
    """
    Process a stream to retrieve noise information
    """
    print(f"\tpreprocess {st_in[0].get_id()}")
    st = st_in.copy()
    # Decimate and then resample the data
    if sampling_rate:
        print(f"\tresample {st[0].stats.sampling_rate} -> {sampling_rate}")
        # Max allowed decimation
        if (st[0].stats.sampling_rate // sampling_rate) > 16:
            st.decimate(16)
        st.resample(sampling_rate)
    st.detrend("demean")
    st.taper(max_percentage=0.05)
    if response:
        try:
            st.remove_response(output="DISP", inventory=inv)
        except ValueError:
            return None
    st.filter("bandpass", freqmin=1/500, freqmax=1/50)
    st.detrend("demean")
    st.taper(max_percentage=0.05)

    return st


def plot_waveforms(st, show=True, save="", tanno=None):
    """
    Plot stream waveforms
    """
    f, ax = plt.subplots()

    for i, tr in enumerate(st):
        line, = plt.plot(tr.times(), tr.data) 
        if tanno:
            plt.text(tanno, tr.data[np.where(tr.times()==tanno)], s=tr.id)
    plt.xlabel("Time (s)")
    plt.ylabel("Displacement (m)")
    plt.grid()
    if show:
        plt.show()
    if save:
        plt.savefig(save)
    plt.close()


def gather_data(choice, comp, sampling_rate=2, response=True): 
    """
    Wrapping the data gathering step into a function
    """
    # Choose time frame
    end = None
    if choice == "noise":
        start = UTCDateTime("2018-055T00:00:00")
        end = start + (24 * 60 * 60) - 1
    elif choice == "alaska":
        start = UTCDateTime("2018-01-23T09:32:00")
    elif choice == "chiapas":
        start = UTCDateTime("2017-09-08T04:49:46.0")
    elif choice == "png":
        start = UTCDateTime("2018-02-25T17:45:08.6")
    elif choice == "mexico":
        start = UTCDateTime("2017-09-19T18:14:48.2")
    elif choice == "local":
        start = UTCDateTime("2018-02-18T07:38:00")
        end = start + (10 * 60)
    
    # Standard end time for teleseismics
    if not end:
        end = start + (3 * 60 * 60)

    # Define full station names
    stations = np.loadtxt("STATIONS", usecols=[1,0], dtype=str)
    stations = [f"{s[0]}.{s[1]}.10.HH{comp}" for s in stations]

    # Stream A
    st_a = Stream()
    for gn in stations:
        try:
            st_a_, inv_a_ = geonet_waveforms(gn, start, end)
        except (FileNotFoundError, ValueError):
            continue
        st_a_ = preprocess(st_a_, inv_a_, sampling_rate, response)
        if st_a_:
            st_a += st_a_

    return st_a


if __name__ == "__main__":
    for choice in ["png"]: # , "mexico", "chiapas", "png"]:
        for comp in ["Z", "N", "E"]:
            st_a = gather_data(choice, comp, response=True)
            plot_waveforms(st_a, show=True, save=f"geonet_{choice}.png", 
                           tanno=984)
