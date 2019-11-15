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

    axes[1].legend(prop={'loc': 'upper left', 'size': 6})
    axes[2].set_xlabel("Time since origin time (sec)")
    plt.sca(axes[0])


def plot_data(geonet_list, st_beacon, title='', normalize=True, show=True, 
              save=''):
    """
    plot the station data

    :param geonet:
    :param st_b:
    :return:
    """
    # plot the beacon station
    f, axes = plt.subplots(3, sharex=True)
    plot_components(axes, st_beacon, color="k", zorder=30, normalize=normalize)

    # plot geonet stations
    for color, st_a in zip(["r", "b", "g"], geonet_list):
        plot_components(axes, st_a, color=color, zorder=20, normalize=normalize)

    starttime = st_beacon[0].stats.starttime
    plt.title(title)
    year = starttime.year

    if save:
        plt.savefig(save, dpi=100)
    if show:
        plt.show()

    plt.close()


def process(starttime, endtime, t0, t1, t_width=4, plot=False, 
            show=True, save=False, normalize=True, **kwargs):
    """
    main processing
    :return:
    """
    comps = ["N", "E", "Z"]

    # Beacon data for each station
    beacon = []
    for i in range(1, 2, 1):
        station = f"RD{i:0>2}"
        st_b = beacon_waveforms(station_name=station, start=starttime,
                                end=endtime, **kwargs)
        # If there is no data, continue
        if not st_b:
            print(f"\tNo data for {station}")
            continue
        else:
            print(f"\t{station}")
        beacon.append(st_b)
    if not beacon:
        return

    # Get GeoNet data for chosen stations
    geonet = []
    for station in ["BKZ", "PXZ", "TSZ"]:
        st_a = geonet_waveforms(station_code=f"NZ.{station}.10.HH?",
                                start=starttime, end=endtime)
        geonet.append(st_a)
        print(f"\t{station}")

    # Preprocess for frequency bands with some width
    for t in range(t0, t1, 1):
        print(f"{t}s")
        if t_width:
            t_min = t0 - t_width / 2
            t_max = t0 + t_width / 2
        else:
            t_min = t0
            t_max = t1

        # Filter geonet data
        geonet_filtered = []
        for st_g in geonet:
            geonet_filtered.append(preprocess(st_g, t_min, t_max))
        
        # Filter beacon data, compare to geonet data
        beacon_filtered = []
        for st_b in beacon:
            st_bf = preprocess(st_b, t_min, t_max)
            bf_netsta = ".".join(st_bf[0].get_id().split('.')[:2])
            print(f"\t{bf_netsta}")
            for st_gf in geonet_filtered:
                gf_netsta = ".".join(st_gf[0].get_id().split('.')[:2])
                print(f"\t\t{gf_netsta}", end=" ")
                for comp in comps:
                    print(f"{comp}:", end=" ")
                    # Max amplitude and time for beacon station
                    tr_bf = st_bf.select(component=comp)[0]
                    amp_max_bf = tr_bf.data.max()
                    time_max_bf = np.where(tr_bf.data == amp_max_bf)[0][0]
                    
                    # Max amplitude and time for geonet station
                    tr_gf = st_gf.select(component=comp)[0]
                    amp_max_gf = tr_gf.data.max()
                    time_max_gf = np.where(tr_gf.data == amp_max_gf)[0][0]

                    # Get amplitude ratio, time difference. Amplitude ratio with
                    # beacon on denom cause we want a scale factor
                    time_diff = time_max_gf - time_max_bf
                    amp_ratio = amp_max_gf / amp_max_bf
                    print(f"A:{amp_ratio:.2E}/T:{time_diff:.2f}s", end=" ")
                print("")
            
            # Plot the data if necessary
            if plot:
                fid_out = (f"./figures/{starttime.year}_{starttime.julday:0>3}_"
                           f"{bf_netsta}_{t_min}_{t_max}.png")
                title = (f"{starttime} T=[{t_min}, {t_max}]s\n{bf_netsta}")
                plot_data(geonet_list=geonet_filtered, st_beacon=st_bf, 
                          title=title, show=show, save=fid_out,
                          normalize=normalize)
        if not t_width:
            return    


def pathing():
    """
    Return the correct paths
    """
    path = "/scale_wlg_persistent/filesets/project/nesi00263/bchow/seismic/" \
           "mseeds/{year}/XX/{sta}/HH?.D/"
    inv_path = "./beacon.xml"
    if not os.path.exists(inv_path):
        # GNS
        path = "/seis/prj/fwi/bchow/data/mseeds/BEACON/{year}/XX/{sta}/HH?.D/"
        inv_path = "/seis/prj/fwi/bchow/data/mseeds/BEACON/DATALESS/beacon.xml"
        if not os.path.exists(inv_path):
            # VUW
            path = "/Users/chowbr/Documents/subduction/seismic/mseeds/" \
                   "BEACON/{year}/XX/{sta}/HH?.D/"
            inv_path = "/Users/chowbr/Documents/subduction/seismic/mseeds/" \
                       "BEACON/DATALESS/beacon.xml"

    return path, inv_path

    
if __name__ == "__main__":
    path, inv_path = pathing()

    # Hand selected earthquakes that produced teleseismics captured by Beacon
    origin_times = [("2018-02-18T07:43:48.0", "2018p130600"),  # 2018p130600 M5.15
                    ("2017-09-08T04:49:46.0", "chiapas"),  # chiapas, mw8.2
                    ("2018-01-23T09:32:00.0", "alaska"),  # alaska, mw7.9
                    ("2018-02-25T17:45:08.6", "png"),   # png, mw7.5
                    ("2017-09-19T18:14:48.2", "mexico"),  # central mexico, mw7.1
                    ]

    # Define the period range for filtering
    t0, t1 = 10, 30
    t_width = None
    for event in origin_times[1:]:
        start, end = event_information(start=event[0], start_pad=0,
                                       end_pad=60*100)
        print(f"{event[1]}: {event[0]}")
        # Filter by small bandwithds to look at frequency dependence
        process(start, end, path=path, inv_path=inv_path,
                plot=True, t0=t0, t1=t1, t_width=t_width, 
                save=True, show=False)










