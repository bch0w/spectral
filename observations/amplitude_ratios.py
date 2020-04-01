"""
Quick script to look at BEACON data for large teleseismic events to check
amplitudes of isntruments in the network
"""
import os
import sys
import glob
import json
import warnings
import numpy as np
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy.signal.cross_correlation import xcorr_3c
from obspy import UTCDateTime, read, Stream, read_inventory

import matplotlib as mpl
mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = .75
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 3


# ignore warnings
warnings.filterwarnings("ignore")


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

    axes[1].legend(loc='upper left', prop={'size': 6})
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

    plt.close("all")


def process(starttime, endtime, t0, t1, t_width=4, plot=False, 
            show=True, save=False, normalize=True, event_name='', **kwargs):
    """
    main processing
    :return:
    """
    from pyatoa.utils.tools.process import trimstreams

    comps = ["N", "E", "Z"]

    # Beacon data for each station
    beacon = []
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

    # Filter beacon data, compare to geonet data
    beacon_filtered = []
    ratio_dict = {}
    for st_b in beacon:
        b_netsta = ".".join(st_b[0].get_id().split('.')[:2])
        print(f"{b_netsta}")
        ratio_dict[b_netsta] = {}
        ratio_dict["periods"] = []
        
        # Preprocess for frequency bands with some width
        for t in range(t0, t1 + 1, 1):
            print(f"{t}s", end=" ")
            if t_width:
                t_min = t - t_width / 2
                t_max = t + t_width / 2
                ratio_dict["periods"].append(t)
            else:
                t_min = t0
                t_max = t1
            
            # Filter Beacon data
            st_bf = preprocess(st_b, t_min, t_max)

            # Filter geonet data
            geonet_filtered = []
            for st_g in geonet:
                st_gf = preprocess(st_g, t_min, t_max)
                geonet_filtered.append(st_gf)
                g_netsta = ".".join(st_gf[0].get_id().split('.')[:2])
                print(f"{g_netsta}", end=" ")
                
                # Initiate dictionary
                if g_netsta not in ratio_dict[b_netsta]:
                    ratio_dict[b_netsta][g_netsta] = {}

                # Loop through each component of both streams
                for comp in comps:
                    print(comp, end=" ")
                    if comp not in ratio_dict[b_netsta][g_netsta]:
                        ratio_dict[b_netsta][g_netsta][comp] = {"time_diff":[],
                                                                "amp_ratio":[]}
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
                    time_diff = float(time_max_gf - time_max_bf) * \
                                                       tr_bf.stats.sampling_rate
                    amp_ratio = amp_max_gf / amp_max_bf
                    
                    ratio_dict[b_netsta][g_netsta][comp]["time_diff"].append(
                                                                    time_diff)
                    ratio_dict[b_netsta][g_netsta][comp]["amp_ratio"].append(
                                                                    amp_ratio)
            print("")
            # Plot the data if necessary
            if plot:
                fid_out = (f"./figures/{event_name}_{b_netsta}_"
                           f"{int(t_min)}_{int(t_max)}.png")
                title = (f"{starttime} T=[{t_min}, {t_max}]s\n{b_netsta}")
                plot_data(geonet_list=geonet_filtered, st_beacon=st_bf, 
                          title=title, show=show, save=fid_out,
                          normalize=normalize)
        if not t_width:
            return    

    return ratio_dict


def plot_ratio_dict(event, path_to_jsons="./"):
    """
    Plot the ratio dictionaries
    """   
    cmarkers = {"N":"o", "E":"s", "Z":"^"}
    gcolors = {"NZ.BKZ":"coral", "NZ.TSZ":"deepskyblue", 
               "NZ.PXZ":"mediumseagreen"}
   
    norm = mpl.colors.Normalize(vmin=1, vmax=22) 
    cmap = mpl.cm.get_cmap('jet')

    for fid in glob.glob(os.path.join(path_to_jsons, f"*{event}*.json")):
        with open(fid) as f:
            ratio_dict = json.load(f)   
        periods = ratio_dict["periods"]
        for sta_b in ratio_dict:
            if sta_b == "periods":
                continue    
            # manual exclusions
            if "chiapas" in fid and "14" in sta_b:  # spike in data
                continue
            b_num = int(sta_b.split("RD")[1])
            for sta_g in ratio_dict[sta_b]:
                for comp in ratio_dict[sta_b][sta_g]:
                    if comp != "Z":
                        continue
                    # time_diff = ratio_dict["time_diff"]
                    amp_ratio = ratio_dict[sta_b][sta_g][comp]["amp_ratio"]
                    plt.scatter(periods, amp_ratio, c='None',
                                edgecolor=gcolors[sta_g], marker=cmarkers[comp], 
                                zorder=10)

            # Plot a line connecting the markers
            plt.plot(periods, amp_ratio, markersize=0, label=sta_b,
                     zorder=5, c=cmap(norm(b_num), alpha=0.5), linewidth=2)
            
            # Write annotation for station and component
            anno = f"{sta_b} {sta_g} {comp}"
            plt.annotate(xy=(periods[-1], amp_ratio[-1]), s=anno, 
                         fontsize=5)

        # plt.legend() 
        # plt.ylim([0,200])
        plt.xlim([5, 40])
        plt.xlabel("Periods (s)")
        plt.yscale("log")
        plt.ylabel("Amplitude Ratio (Geonet/Beacon)")
        plt.title(f"BKZ (r), TSZ (b), PXZ (g)\n")
        plt.show()


def plot_by_beacon(sta_b="XX.RD01", event='', path_to_jsons="./", 
                   overwrite=False):
    """
    Plot the ratio dictionaries
    """   
    cmarkers = {"N":"o", "E":"s", "Z":"^"}
    gcolors = {"NZ.BKZ":"coral", "NZ.TSZ":"deepskyblue", 
               "NZ.PXZ":"mediumseagreen"}
    ecolors = {"chiapas": "coral", "alaska": "deepskyblue", 
               "png": "mediumseagreen", "mexico": "fuchsia"}
    glinestyles = {"NZ.BKZ":"--", "NZ.TSZ":"-", "NZ.PXZ":"-."}
     
    periods_check = None
    fids = glob.glob(os.path.join(path_to_jsons, f"*{event}*.json"))
    for fid in fids:
        event_name = os.path.basename(fid.split('_')[0])
        print(event_name)
        with open(fid) as f:
            ratio_dict = json.load(f)   
    
        # make sure the periods are the same between files
        periods = ratio_dict["periods"]
        if not periods_check:
            periods_check = periods
        else:
            assert(periods_check == periods)
       
        # make sure the beacon station actually has data 
        if not sta_b in ratio_dict.keys():
            print(f"{sta_b} not in ratio_dict")
            continue
        
        # make sure we haven't run this before
        fig_title = f"./figures/{sta_b}.png"
        # if os.path.exists(fig_title) and not overwrite:
        #     print("figure exists")
        #     continue

        # get amplitude ratios for each geonet station listed under this station
        for sta_g in ratio_dict[sta_b].keys():
            for comp in ratio_dict[sta_b][sta_g]:
                amp_ratio = ratio_dict[sta_b][sta_g][comp]["amp_ratio"]
                
                # plot color by station/component
                plt.scatter(periods, amp_ratio, c='None',
                            edgecolor=ecolors[event_name], marker=cmarkers[comp], 
                            zorder=10)
                plt.plot(periods, amp_ratio, markersize=0, zorder=5, 
                         c=ecolors[event_name], linewidth=1, alpha=0.5,
                         linestyle=glinestyles[sta_g])

        # Plot a line connecting the marker
        plt.plot(periods, amp_ratio, markersize=0, zorder=5, 
                 c=ecolors[event_name], linewidth=2, 
                 linestyle=glinestyles[sta_g], 
                 label=f"{event_name} {sta_g} {comp}")

                # # plot color by station/component
                # plt.scatter(periods, amp_ratio, c='None',
                #             edgecolor=gcolors[sta_g], marker=cmarkers[comp], 
                #             zorder=10)
                # # Plot a line connecting the marker
                # if comp == "Z":
                #     plt.plot(periods, amp_ratio, markersize=0, zorder=5, 
                #              c=gcolors[sta_g], linewidth=2, label=f"{sta_g}_Z")
                # else:
                #     plt.plot(periods, amp_ratio, markersize=0, zorder=5, 
                #              c=gcolors[sta_g], linewidth=1, alpha=0.5,
                #              linestyle="--")


    plt.legend(loc="upper left")
    plt.xlabel("Periods (s)")
    # plt.yscale("log")
    plt.grid(True, linestyle="--")
    plt.ylabel("Amplitude Ratio (Geonet/Beacon)")
    plt.title(f"{sta_b}\n"
              "BKZ --, TSZ -, PXZ -.; "
              "North (o), East (square), Z (triangle)")

    plt.savefig(fig_title, figsize=(8,8))
    # plt.show()
    plt.close("all")

 
if __name__ == "__main__":
    process_events = True
    path, inv_path = pathing()

    # Hand selected earthquakes that produced teleseismics captured by Beacon
    origin_times = [("2018-02-18T07:43:48.0", "2018p130600"),  # ml5.2
                    ("2017-09-08T04:49:46.0", "chiapas"),      # mw8.2
                    ("2018-01-23T09:32:00.0", "alaska"),       # mw7.9
                    ("2018-02-25T17:45:08.6", "png"),          # mw7.5
                    ("2017-09-19T18:14:48.2", "mexico"),       # mw7.1
                    ]
   
    # processing saved into jsons so it should only be run once 
    if process_events:
        for event in origin_times:
            t0, t1 = 1, 6
            t_width = 4
            start, end = event_information(start=event[0], start_pad=0,
                                           end_pad=60*100)
            json_fid = f"{event[1]}_{t0}_{t1}_{t_width}.json"
            if not os.path.exists(json_fid):
                print(f"{event[1]}: {event[0]}")
                # Filter by small bandwithds to look at frequency dependence
                ratio_dict = process(start, end, path=path, inv_path=inv_path,
                                     plot=True, t0=t0, t1=t1, t_width=t_width, 
                                     save=True, show=False, event_name=event[1])
                with open(json_fid, "w") as f:
                    json.dump(ratio_dict, f, indent=4)
            else:
                with open(json_fid, "r") as f:
                    ratio_dict = json.load(f)

    # plot the data
    for i in range(1, 24, 1):
        station = f"XX.RD{i:0>2}"
        print(station)
        plot_by_beacon(station, path_to_jsons="./jsons")

    # for event in origin_times:
    #     plot_ratio_dict(event=event[1], path_to_jsons="./") 








 
