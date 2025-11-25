"""
Implements the Segmented Window detector as described in Tian et al. (2022)
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime

def set_plot_aesthetic(
        ax, ytick_fontsize=9., xtick_fontsize=9., tick_linewidth=1.5,
        tick_length=5., tick_direction="in", ytick_format="sci",
        xlabel_fontsize=10., ylabel_fontsize=10., axis_linewidth=1.5,
        spine_zorder=8, title_fontsize=10.,
        spine_top=True, spine_bot=True, spine_left=True, spine_right=True,
        xtick_minor=None, xtick_major=None, ytick_minor=None, ytick_major=None,
        xgrid_major=True, xgrid_minor=True, ygrid_major=True, ygrid_minor=True,
        **kwargs):
    """
    Set a uniform look for figures, stolen from PySEP
    """
    ax.title.set_fontsize(title_fontsize)
    ax.tick_params(axis="both", which="both", width=tick_linewidth,
                        direction=tick_direction, length=tick_length)
    ax.tick_params(axis="x", labelsize=xtick_fontsize)
    ax.tick_params(axis="y", labelsize=ytick_fontsize)
    ax.xaxis.label.set_size(xlabel_fontsize)
    ax.yaxis.label.set_size(ylabel_fontsize)

    # Thicken up the bounding axis lines
    for axis, flag in zip(["top", "bottom", "left", "right"],
                          [spine_top, spine_bot, spine_left, spine_right]):
        # Deal with the case where command line users are inputting strings
        if isinstance(flag, str):
            flag = bool(flag.capitalize() == "True")
        ax.spines[axis].set_visible(flag)
        ax.spines[axis].set_linewidth(axis_linewidth)

    # Set spines above azimuth bins
    for spine in ax.spines.values():
        spine.set_zorder(spine_zorder)

    # Scientific format for Y-axis
    if ytick_format == "sci":
        try:
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        except AttributeError:
            # If we are in log format axis this will not work
            pass
    else:
        try:
            ax.ticklabel_format(axis="y", style=ytick_format)
        except AttributeError:
            pass

    # Set xtick label major and minor which is assumed to be a time series
    if xtick_major:
        ax.xaxis.set_major_locator(MultipleLocator(float(xtick_major)))
    if xtick_minor:
        ax.xaxis.set_minor_locator(MultipleLocator(float(xtick_minor)))
    if ytick_minor:
        ax.yaxis.set_major_locator(MultipleLocator(float(ytick_major)))
    if ytick_major:
        ax.yaxis.set_minor_locator(MultipleLocator(float(ytick_minor)))

    plt.sca(ax)
    if xgrid_major:
        plt.grid(visible=True, which="major", axis="x", alpha=0.2, linewidth=1)
    if xgrid_minor:
        plt.grid(visible=True, which="minor", axis="x", alpha=0.2, linewidth=.5)
    if ygrid_major:
        plt.grid(visible=True, which="major", axis="y", alpha=0.2, linewidth=1)
    if ygrid_minor:
        plt.grid(visible=True, which="minor", axis="y", alpha=0.2, linewidth=.5)


def segmented_window(data, window_length_s, sampling_rate):
    """
    Segmented window detector, provides a ratio time series and sample indices
    for detections within the `data` 

    :type data: np.array
    :param data: seismic signal, arbitrary amplitudes but counts preferred to
        mimic that we won't be able to remove response from our instrument
    :type window_length: int
    :param window_length: length in seconds of time window used for defining 
        ratio
    :type sampling_rate: float
    :param sampling_rate: sampling rate of the data
    :rtype:
    :return:
    """
    # Set some needed values
    window_length_samples = int(window_length_s * sampling_rate)
    nsamples = len(data)
    total_windows_allowed = nsamples // window_length_samples
    
    # Loop through data stream one window at a time
    #   i: current window starttime index
    #   j: current window endtime index
    #   n: current window index
    #   m: previous window
    i = 0
    ratios = []
    for n, i in enumerate(np.arange(0, nsamples, window_length_samples, 
                                    dtype=int)):
        # We will have unfinished windows, append zeros to fill
        if n >= total_windows_allowed:
            ratios += [0] * (nsamples - len(ratios))
            break

        j = i + window_length_samples 

        # Calculate the window average which will be used as the previous window
        win_m_avg = np.sum(np.abs(data[i:j])) / window_length_samples

        # If this is the first window, move on
        if n == 0:
            ratios[i:j] = np.zeros(window_length_samples)
            continue

        # For all windows n > 0, start calculating the ratio
        for sample in np.arange(i, j, 1):
            ratio = np.abs(data[sample]) / win_m_avg
            ratios.append(ratio)

    return ratios

def count_threshold(ratios, R, val=1):
    """
    Return a time series that is 0 or `val` if a threshold is past
    """
    thresholds = []
    counts = 0
    for r in ratios:
        if r >= R:
            thresholds.append(val)
            counts += 1
        else:
            thresholds.append(-0.1)

    return thresholds, counts


if __name__ == "__main__":
    window_length_s = 120
    R = 7
    if False:
        freqmin=2  # 0.01
        freqmax=8  # 10
    else:
        freqmin=.1
        freqmax=10

    lw = 1

    stations = [
            f"IU.POHA.00.BHZ.2025.087",
            f"VH.DESD.EX.EHZ.2025.087",
            f"VH.MITD.EX.EHZ.2025.087",
            f"VH.KAED.EX.EHZ.2025.087",
            ]



    stations = [
            f"VH.KAED.EX.EHZ.2025.087",
            f"VH.KAED.CV.EHZ.2025.087",
            f"VH.KAED.WS.EHZ.2025.087",
            f"VH.DESD.EX.EHZ.2025.087",
            f"VH.DESD.CV.EHZ.2025.087",
            f"VH.DESD.WS.EHZ.2025.087",
            ]

    f, axs = plt.subplots(3, sharex=True, figsize=(8,5), dpi=200)
    for ax in axs:
        set_plot_aesthetic(ax, ytick_format="plain", xgrid_major=False,
                           xgrid_minor=False, ygrid_major=False, 
                           ygrid_minor=False)
    f.subplots_adjust(hspace=0)

    # Read and prep data
    path = "/Users/chow/Work/research/venusseis/HAVO_data_analysis/VEXAG_slide"

    # Actual time range to be processing
    start = UTCDateTime("2025-03-28T06:20:00")
    end = UTCDateTime("2025-03-28T07:00")
    buffer = 60 * 10  # sec

    # TauP expected arrivals
    expected_p = UTCDateTime("2025-03-28T06:34:31") - start
    expected_s = UTCDateTime("2025-03-28T06:46:00") - start

    max_amp, max_rat = 0, 0

    # Calculate offset values
    _buf = 0.025
    _min = 1 - _buf * len(stations)
    vals = np.arange(_min, 1, _buf)
    for offset, station in enumerate(stations):
        net, sta, loc, cha, year, jday = station.split(".")
        st = read(os.path.join(path, station))

        # Trim big
        st.trim(start - buffer, end + buffer)

        # Preproc
        st.taper(0.05)
        st.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, 
                  zerophase=True)

        # Trim tails
        st.trim(start, end)

        # Normalize Amplitudes
        st[0].data /= st[0].data.max()

        # Feed into trigger algo
        ratios = segmented_window(data=st[0].data, 
                                  window_length_s=window_length_s, 
                                  sampling_rate=st[0].stats.sampling_rate)
        # Figure out when we reach "detection"
        thresholds, counts = count_threshold(ratios, R, vals[offset])

        # Plot all the different metrics
        axs[0].plot(st[0].times(), st[0].data + offset * 2, lw=lw, 
                    label=f"{net}.{sta}.{loc}.{cha}")
        axs[1].plot(st[0].times(), ratios, lw=lw)
        axs[2].plot(st[0].times(), thresholds, f"C{offset}o-", lw=0.5,
                    markersize=1, alpha=0.5, zorder=5 + 2*(1/(offset+1)),
                    label=counts)

        if max(ratios) > max_rat:
            max_rat = max(ratios)
        if st[0].data.max() > max_amp:
            max_amp = st[0].data.max()

    # Demarcate the time windows
    tmax = st[0].times().max()
    kwargs = dict(alpha=0.25, ls=":", lw=0.6, c="k")
    for i in np.arange(0, tmax, window_length_s):
        axs[0].axvline(i, **kwargs)
        if i == 0:
            axs[1].axvline(i, label=f"DT={window_length_s}s", **kwargs)
        else:
            axs[1].axvline(i, **kwargs)

    # Plot the threshold
    axs[1].axhline(R, c="k", ls="--", lw=1, label=f"Threshold={R}")

    # Expected Body Wave Arrivals
    axs[0].axvline(expected_p, c="magenta", ls="--", lw=.75, 
                   label="Expected P, S")
    axs[0].axvline(expected_s, c="magenta", ls="--", lw=.75)

    # Set ylims
    # axs[0].set_ylim([-1 * max_amp, max_amp])
    axs[1].set_ylim([0, max_rat + 1])
    axs[2].set_ylim([0, 1.05])

    # Turn off axes I don't need
    axs[0].set_yticks([])
    axs[2].set_yticks([])
    

    plt.xlabel("Time since earthquake [s]")
    plt.xlim([0, tmax])
    axs[0].set_ylabel("Norm. Amp.")
    axs[1].set_ylabel(f"Ratio")
    axs[2].set_ylabel("Detection")
    axs[0].set_title(f"Myanmar Mw7.7 "
                     f"2025-03-28T06:20:52; "
                     f"$\Delta$=99$^\circ$; [{freqmin}, {freqmax}]Hz")

    axs[0].legend(loc="upper right", fontsize=5)
    axs[1].legend(loc="upper right", fontsize=5)
    axs[2].legend(loc="upper right", fontsize=5)
    plt.show()


    # # Store the ratio as an ObsPy stream for easier use later
    # st_ratio = st.copy()
    # st_ratio[0].data = np.array(ratios)
    # st_ratio[0].stats.location = "RT"

