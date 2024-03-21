"""
Look at the amplitude levels for given stations
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

        c = Client("GEONET")
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


def beacon_waveforms(station_name, start, end):
    """
    Get beacon station waveforms based on station number
    """
    print(f"Beacon waveforms {station_name}")
    code = f"{station_name}.D.{start.year}.{start.julday:0>3}"
    net, sta, loc, cha, d, year, jday = code.split(".")

    path = (f"/scale_wlg_persistent/filesets/project/nesi00263/bchow/seismic/" \
            f"mseeds/{year}/XX/{sta}/{cha}.D/")
    inv_path = os.path.expanduser(os.path.join("~", "primer", "auxiliary",
                                               "stationxml", "beacon.xml")
                                  )
    fid = os.path.join(path, code)
    if os.path.exists(fid):
        st_beacon = read(fid)
    else:
        raise FileNotFoundError(f"Beacon data not found for {station_name}"
                                f"at time: {start.year}.{start.julday}\n"
                                f"{path}")

    st_beacon = st_beacon.trim(start, end)
    if len(st_beacon) > 1:
        raise DataGapError()
    inv = read_inventory(inv_path)

    # Will only attach the relevant response
    st_beacon.attach_response(inv)
    
    return st_beacon


def preprocess(st_in, inv=None, sampling_rate=None, response=True, 
               normalize=True):
    """
    Process a stream to retrieve noise information
    """
    tmin, tmax = 10, 30
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
        st.remove_response(output="DISP", inventory=inv)
    st.filter("bandpass", freqmin=1/tmax, freqmax=1/tmin, zerophase=True)
    
    # Special preprocessing for beacon data
    if st[0].stats.network == "XX":
        scale = 75.
        code = int(st[0].stats.station[-2:])
        # 60s instruments
        if code in [1, 2, 4, 5, 6, 7, 8, 9, 13, 14, 15, 18]:
            print(f"scaling {st[0].id} by {scale}")
            for tr in st:
                tr.data *= scale
                tr.stats.location = "60"
        else:
            for tr in st:
                tr.stats.location = "30"

    st.detrend("demean")
    st.taper(max_percentage=0.05)
    
    if normalize:
        for tr in st:
            tr.data /= tr.max()
    t = st[0].stats.starttime
    st.write(f"{st[0].id}.{tmin}-{tmax}.{t.year}-{t.julday}",
             format='mseed')

    return st


def plot_waveforms(st, show=True, save=""):
    """
    Plot stream waveforms
    """
    for i, tr in enumerate(st):
        plt.plot(tr.times(), tr.data, label=tr.get_id())
    plt.xlabel("Time (s)")
    plt.ylabel("Displacement (m)")
    plt.title("Waveform comparisons")
    plt.grid()
    plt.legend()
    if show:
        plt.show()
    if save:
        plt.savefig(save)
    plt.close()


def gather_data(choice, min_freq=1/100, sampling_rate=2, response=True, 
                comp="Z"): 
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
    bcn_fmt = "{:0>2}"  # sneak a format into a formatted string
    geonet = [f"NZ.BFZ.10.HH{comp}", f"NZ.PXZ.10.HH{comp}", 
              f"NZ.TSZ.10.HH{comp}"]
    beacon = f"XX.RD{bcn_fmt}.10.HH{comp}"

    # Stream B for Beacon
    st_b = Stream()
    for i in []:
    # for i in [1, 2, 4, 5, 6, 7, 8, 9, 13, 14, 15, 18]:
        try:
            st_b_ = beacon_waveforms(beacon.format(i), start, end)
        except FileNotFoundError:
            print("No data")
            continue
        except DataGapError:
            print("Non continuous data")
            continue
        st_b_ = preprocess(st_b_, None, sampling_rate, response=response)
        st_b += st_b_
    
    # Stream A for GeoNet
    st_a = Stream()
    for gn in geonet:
        st_a_, inv_a_ = geonet_waveforms(gn, start, end)
        st_a_ = preprocess(st_a_, inv_a_, sampling_rate, response=response)
        st_a += st_a_

    return st_a, st_b


if __name__ == "__main__":
    choices = ["alaska"] 
    choices = ["mexico", "chiapas", "png"]
    for choice in ["alaska"]:
        for comp in ["N"]:  #, "N", "E"]:
            st_a, st_b = gather_data(choice, response=True, comp=comp)
