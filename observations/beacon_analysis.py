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


def preprocess(st_in, freqmin, inv=None, sampling_rate=None, response=True):
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
        st.remove_response(output="DISP", inventory=inv)
    st.filter("highpass", freq=freqmin)
    st.detrend("demean")
    st.taper(max_percentage=0.05)

    return st


def amplitude_ratios(st_beacon, st_geonet, t_min, t_max, t_width=4):
    """
    Filter streams for narrow bandpasses and compare amplitude values
    """
    amp_dict = {}
    central_t = []
    for period_min in range(t_min, t_max + 1, 1):
        st_b = st_beacon.copy()
        st_g = st_geonet.copy()
        period_max = period_min + t_width
        central_t.append((period_max + period_min) / 2)

        # Filter for this narrow bandpass
        st_b.filter("bandpass", freqmin=1 / period_max, freqmax=1 / period_min)
        st_g.filter("bandpass", freqmin=1 / period_max, freqmax=1 / period_min)

        # Compare a few values
        rms_values, max_values = [], []
        for tr in (st_g + st_b):
            if tr.get_id() not in amp_dict.keys():
                amp_dict[tr.get_id()] = {"max":[], "mean":[], "rms":[]}
            amp_dict[tr.get_id()]["max"].append(tr.data.max())
            amp_dict[tr.get_id()]["mean"].append(tr.data.mean())
            amp_dict[tr.get_id()]["rms"].append(np.sqrt(np.mean(tr.data ** 2)))
  

    # Turn the lists into numpy arrays
    for key1, item1 in amp_dict.items():
        for key2, item2 in item1.items():
            item1[key2] = np.array(item2)
    
    amp_dict["periods"] = np.array(central_t) 

    return amp_dict


def plot_waveforms(st, show=True, save=""):
    """
    Plot stream waveforms
    """
    colors = ["k", "r"]
    for i, tr in enumerate(st):
        plt.plot(tr.times(), tr.data, c=colors[i],
                 label=tr.get_id(), alpha=1 / (i + 1))
    plt.xlabel("Time (s)")
    plt.ylabel("Velocity (m/s)")
    plt.title("Waveform comparisons, normalized\n 1/100 to 2Hz SR")
    plt.grid()
    plt.legend()
    if show:
        plt.show()
    if save:
        plt.savefig(save)
    plt.close()


def plot_periodogram(st, freqmin=None, show=True, save=""):
    """
    Create periodogram plots for a given stream
    """
    print("plotting")
    colors = ["k", "r", "b", "orange"]
    compcol = ["gray", "orange"]
    assert (len(colors) > len(st))

    ybounds = []
    bins = []
    for i, tr in enumerate(st):
        print(f"\tperiodogram {tr.get_id()}")
        f, pxx = signal.periodogram(tr.data, tr.stats.sampling_rate)
        plt.plot(f, oxx, c=colors[i], label=tr.get_id(), zorder=50,
                 alpha=1 / (i + 1))

        # Bin data and take means
        bin_x, bin_y = [], []
        for freq_low in np.arange(1 / 200, 1, .005):
            freq_high = freq_low + 0.05
            mean_val = np.mean(Pxx[np.where((f > freq_low) & (f < freq_high))])
            bin_y.append(mean_val)
            bin_x.append(np.mean([freq_low, freq_high]))

        # Plot the mean line through
        plt.plot(bin_x, bin_y, marker='o', markeredgecolor='w',
                 color=compcol[i], zorder=100, alpha=1 / (i + 1))
        bins.append(bin_y)

        # Retrieve min and max values for bounds
        ybounds.append(pxx[np.where(f > freqmin)].min())
        ybounds.append(pxx[np.where(f > freqmin)].max())

    # Plot dominant periods for reference
    for freqs in [freqmin, 1 / 50, 1 / 20, 1 / 10, 1 / 5, 1 / 2]:
        plt.axvline(x=freqs, linestyle='--', zorder=25, c="k", alpha=0.5)
        plt.text(freqs, min(ybounds), f"{int(1/freqs):d}")

    # Only bound the data that were interested in
    plt.xlim([freqmin or 0, 1 / 2])
    plt.ylim([min(ybounds), max(ybounds)])

    # Log scale for better look at period ranges
    plt.yscale("log")
    plt.xscale("log")

    plt.xlabel("Frequency (Hz)")
    plt.ylabel("PSD")
    start = st[0].stats.starttime
    end = st[0].stats.endtime

    # Tell user amplitude ratios
    print("amplitude ratios")
    for i in range(len(bin_x)):
        print(f"{bin_x[i]}: {bins[1][i] / bins[0][i]:.2f}")

    plt.title(f"Periodogram Comparisons "
              f"{freqmin} to {st[0].stats.sampling_rate} Hz\n"
              f"{start.year}.{start.julday:0>3} - {end.year}.{end.julday:0>3}\n"
              f"{(end-start)/(60*60):.2f} hours data")
    plt.legend()
    if save:
        plt.savefig(save)
    if show:
        plt.show()
    plt.close()


def plot_amplitude_ratio(st_a, st_b, amp_dict, show=True, save=True, 
                         direct=False, **kwargs):
    """
    Plot the ratio of ampltiudes
    """
    choice = kwargs.get("choice", "")    
    style = kwargs.get("style", "")

    def plot_attributes():
        plt.xlabel("Period (s)")    
        plt.legend(loc="center", ncol=4, fontsize="x-small", framealpha=0.25)
        plt.yscale("log")
        plt.grid()

    # Set the colorbar
    vmax = len(st_a + st_b)
    norm = mpl.colors.Normalize(vmin=1, vmax=vmax)
    cmap = mpl.cm.get_cmap("tab20b")

    # Plot direct amplitdues for given style, e.g. rms
    if direct:
        for i, tr in enumerate(st_a + st_b):
            # Overwrite for NZ statinos, ensure consistent coloring
            if "NZ" in tr.get_id():
                color = "k"
                linestyle= '--'
            else:
                y = int(tr.get_id().split('.')[1][2:])
                color = cmap(norm(y))
                linestyle='-'
                
            plt.plot(amp_dict["periods"], amp_dict[tr.get_id()][style], 
                     label=tr.get_id(), c=color, linestyle=linestyle)
        plt.title(f"Amplitudes {choice} {style}\n"
                  f"{tr.stats.starttime} - {tr.stats.endtime}")
        plt.ylabel("RMS Displacement (m)")
        plot_attributes()
        if save:
            plt.savefig(f"./figures/"
                        f"amplitudes_{choice}_{style}.png")
        if show:
            plt.show()   
    else: 
        geonet = []
        for tr in st_a:
            geonet.append(tr.get_id())

        for g in range(len(geonet)):
            baseline = amp_dict[geonet[g]][style]
            for i, tr in enumerate(st_a + st_b):
                ratio = baseline / amp_dict[tr.get_id()][style]
                #  Check the marker types based on the ratio
                if ratio.mean() <= 50:
                    linestyle = 'x--'
                else:
                    linestyle = 'o-'
                # Ensure consistent coloring
                if "NZ" in tr.get_id():
                    y = vmax
                else:
                    y = int(tr.get_id().split('.')[1][2:])
                plt.plot(amp_dict["periods"], ratio, linestyle,  markersize=3, 
                         label=tr.get_id(), c=cmap(norm(i)))
            plt.title(f"Amplitude {style} w.r.t {geonet[g]}\n"
                      f"{tr.stats.starttime} - {tr.stats.endtime}")
            plt.ylabel("RMS Amplitude Ratios")
            plot_attributes()
            if save:
                plt.savefig(f"./figures/"
                            f"ratios_{geonet[g]}_{choice}_{style}.png")
            if show:
                plt.show()   


def gather_data(choice, min_freq=1/100, sampling_rate=2, response=True): 
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
    geonet = ["NZ.BFZ.10.HHZ", "NZ.PXZ.10.HHZ", "NZ.TSZ.10.HHZ"]
    beacon = "XX.RD{:0>2}.10.HHZ"

    # Stream B
    st_b = Stream()
    for i in range(1, 23):
        try:
            st_b_ = beacon_waveforms(beacon.format(i), start, end)
        except FileNotFoundError:
            print("No data")
            continue
        except DataGapError:
            print("Non continuous data")
            continue
        st_b_ = preprocess(st_b_, min_freq, None, sampling_rate, response)
        st_b += st_b_

    # Stream A
    st_a = Stream()
    for gn in geonet:
        st_a_, inv_a_ = geonet_waveforms(gn, start, end)
        st_a_ = preprocess(st_a_, min_freq, inv_a_, sampling_rate, response)
        st_a += st_a_

    amp_dict = amplitude_ratios(st_a, st_b, t_min=1, t_max=60)

    return st_a, st_b, amp_dict


if __name__ == "__main__":
    choices = ["alaska", "chiapas", "png", "mexico"]
    teleseismics = ["alaska", "mexico", "chiapas", "png"]
    amp_dicts = []
    for choice in choices:
        st_a, st_b, amp_dict = gather_data(choice, response=True)
        amp_dicts.append(amp_dict) 
        # plot_amplitude_ratio(st_a, st_b, amp_dict, style="rms", choice=choice,
        #                      show=True, save=True, direct=False)

    # for key in amp_dicts[0].keys():
    #     if key == "periods":
    #         continue
    #     i = 0
    #     ratios = np.zeros(len(amp_dicts[0]["periods"]))
    #     for event in amp_dicts:
    #         try:
    #             station = event[key]['rms']
    #             baseline = event["NZ.BFZ.10.HHZ"]['rms']
    #             ratio = baseline / station
    #             ratios += (baseline / station)
    #             i += 1 
    #         except KeyError:
    #              continue
    #     threesec = np.where(amp_dicts[0]["periods"] == 3.)[0]
    #     sixsec = np.where(amp_dicts[0]["periods"] == 6.)[0]
    #     tensec = np.where(amp_dicts[0]["periods"] == 10.)[0]
    #     
    #     print(f"{key}: 3s={ratios[threesec] / i}, 6s={ratios[sixsec] / i}, "
    #           f"10s={ratios[tensec] / i}")
    #     # plt.plot(amp_dicts[0]["periods"], ratios / i)
    #     # plt.show()

    # import ipdb;ipdb.set_trace();

    # # Create plots
    # fid_out = os.path.join("figures",
    #                        f"{geonet.split('.')[1]}.{beacon.split('.')[1]}."
    #                        f"{start.year}.{start.julday:0>3}"
    #                        "_{}.png"
    # )
    # plot_periodogram(st_a + st_b, min_freq, show=False,
    #                  save=fid_out.format("periodogram"))
    plot_waveforms(st_a + st_b, show=False, save=fid_out.format("waveform"))
