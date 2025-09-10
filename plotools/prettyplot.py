"""
Pretty plotting of ObsPy Streams, can either be ObsPy 'read'able format or
two-column ASCII output for SPECFEM synthetics. Most things controlled by parser

.. rubric

    pplot <STA> --xlim 0 200 --fmin 1 --fmax 4 --taper 0.5 --c k --lw 0.5

    pplot NK6_IC.MDJ.CXZ.sem.ascii --differentiate 1 --ylabel 'Velocity [m/s]' \
            -t0 50 --fmin 2 --fmax 4 --save NK6.png --xlim 0 200 \
            --taup Pn Pg Sn Sg --dist 378 --depth 0

    pplot EQ2_IC.MDJ.CXZ.sem.ascii --differentiate 1 --ylabel 'Velocity [m/s]' \
            -t0 50 --fmin 2 --fmax 4 --save EQ2.png --xlim 0 200 \
            --taup Pn Pg Sn Sg --dist 378 --depth 3 --ylim -.0004 .0004 \
            -ta 'SPECFEM3D_GLOBE s20rts_crust1'
"""
import argparse
import sys
import os
import math
import matplotlib.pyplot as plt
import numpy as np

from dateutil.rrule import MINUTELY, SECONDLY
from matplotlib import mlab
from matplotlib.colors import Normalize
from matplotlib.dates import date2num, AutoDateLocator
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

from obspy import read, UTCDateTime, Stream
from obspy.taup import TauPyModel
from obspy.geodetics import kilometers2degrees
from obspy.imaging.cm import obspy_sequential
from obspy.imaging.util import ObsPyAutoDateFormatter

from pysep import read_sem


SECONDS_PER_DAY = 3600.0 * 24.0


def parse_args():
    """All modifications are accomplished with command line arguments"""
    parser = argparse.ArgumentParser()

    # Waveform Processing
    parser.add_argument("fid", nargs="+", help="required, file ID(s)")
    parser.add_argument("-tp", "--taper", nargs="?", type=float, default=0,
                        help="optional taper percentange")
    parser.add_argument("-f1", "--fmin", nargs="?", type=float, default=None,
                        help="optional filtering freqmin in Hz")
    parser.add_argument("-f2", "--fmax", nargs="?", type=float, default=None,
                        help="optional filtering freqmax in Hz")
    parser.add_argument("-z", "--zerophase", action="store_true", default=False,
                        help="apply zerophase filter or not")
    parser.add_argument("-C", "--corners", type=int, default=4,
                        help="number of corners to apply on ")
    parser.add_argument("-r", "--resample", type=float, default=None,
                        help="resample the data prior to processing")
    parser.add_argument("-t0", "--t0", nargs="?", type=float, default=0,
                        help="SPECFEM USER_T0 if synthetics")
    parser.add_argument("-ts", "--tstart", nargs="?", type=float, default=0,
                        help="for relative time axis, set t0 value, defaults 0")
    parser.add_argument("--integrate", nargs="?", type=int, default=0,
                        help="Integrate the time series, will demean and taper,"
                             " value for integrate will be the number of times")
    parser.add_argument("--differentiate", nargs="?", type=int, default=0,
                        help="Diff. the time series, will demean and taper,"
                             " value for integrate will be the number of times")
    parser.add_argument("--trim_pct", nargs="?", type=float, default=1.0,
                        help="Buffer amount in time to trim the data if using " 
                             "`xlim` to prevent processing the tails of the " 
                             "cut data, given in percentage. That is if the " 
                             "requested cut length is 30s, an additional 30s " 
                             "on either side of the trace will be cut and " 
                             "processed")

    # TauP Phase Arrivals
    parser.add_argument("--tp_phases", nargs="+", type=str, default=None,
                        help="taup phase arrivals, requires all 'tp_*' params")
    parser.add_argument("--tp_model", nargs="?", type=str, default="iasp91",
                        help="taup model, defaults to 'iasp91'")
    parser.add_argument("--tp_dist", nargs="?", type=float, default=None,
                        help="TauP source receiver distance in km")
    parser.add_argument("--tp_depth", nargs="?", type=float, default=None,
                        help="TauP source depth km")
    
    # Spectrogram (all parameters end with _s)
    parser.add_argument("--spectrogram", action="store_true", default=False,
                        help="plot spectrogram of the raw trace. See all " 
                             "'sp_*' parameters to control the look")
    parser.add_argument("--sp_cmap", nargs="?", type=str, 
                        default="nipy_spectral",
                        help="colormap of the spectrogram")
    parser.add_argument("--sp_numcol", nargs="?", type=int, default=256,
                        help="number of colors in colormap of the spectrogram")
    parser.add_argument("--sp_dbscale", action="store_true",
                        help="use dB scaling for colors/amplitudes in spectro")    
    parser.add_argument("--sp_logscale", action="store_true",
                        help="turn on log scale for spectrogram y-axis")
    parser.add_argument("--sp_idx", type=int, default=0,
                        help="Iff multiple waveforms plotted, choose which of "
                             "them to use for the spectrogram. Defaults to 0. " 
                             "order is alphabetical")
    
    # Stream Gauge (VERY CUSTOM, ONLY FOR GULKANA EXPERIEMENT)
    parser.add_argument("--stream_gage", action="store_true", default=False,
                        help="For GULKANASEIS data only, plots stream gauge " 
                             "data at the bottom of the waveform plot with a " 
                             "twin X axis")
    parser.add_argument("--sg_rel", action="store_true", default=False,
                        help="For GULKANASEIS data only, make y-axis relative "
                            "so that the lowest plotted value is 0.")
    parser.add_argument("--sg_units", type=str, default="m",
                        help="For GULKANASEIS data only, choose units to plot " 
                             "stream gage data. Options: ft, in, m, cm")

    # Plot Aesthetics
    parser.add_argument("-x", "--xlim", nargs="+", default=None,
                        help="time axis limits in s or if `time`=='a' then "
                             "values should be in datetime, see tmarks")
    parser.add_argument("-y", "--ylim", nargs="+", type=float, default=None,
                        help="amplitude axis limits in s")
    parser.add_argument("-t", "--time", nargs="?", type=str, default="s",
                        help="units for x-axis/time axis. choice: 's'econds "
                             "(default), 'm'inutes, 'h'ours, 'a'bsolute (wip)."
                             "If using 'a' you may add '+i' or '-i' to "
                             "time shift the array, e.g., to go from UTC to "
                             "local time. E.g., 'a-7' will subtract 7 hours.")
    parser.add_argument("-c", "--colors", nargs="+", type=str, default="k",
                        help="color of the time series line, number of inputs "
                             "must match the length of `fid`")
    parser.add_argument("-l", "--labels", nargs="+", type=str, default=None,
                        help="optional labels legend, must match len of `fid`")
    parser.add_argument("-lw", "--linewidth", nargs="?", type=float, 
                        default=0.5, help="linewidth of the time series line")
    parser.add_argument("--ylabel", nargs="?", type=str, default="Vel. [m/s]",
                        help="label for units, defaults to displacement")
    parser.add_argument("--title", nargs="?", type=str, default=None,
                        help="title of the figure, defaults to ID and fmin/max")
    parser.add_argument("-ta", "--title_append", nargs="?", type=str, 
                        default="", help="append to default title")
    parser.add_argument("--maxticks", type=int, default=6, 
                        help="max ticks if --time='a'")
    
    # Time Marks
    parser.add_argument("-tm", "--tmarks", nargs="+", 
                        help="plot vertical lines at given relative times, "
                             "should match the units of `time`. If `time`=='a' "
                             "then each tmark should be a datetime "
                             "YYYY-MM-DDTHH:MM:SS")
    parser.add_argument("--tmarks_c", nargs="+", default="k",
                        help="colors for each of the time marks, should "
                             "either be single letter for all marks or match " 
                             "length of tmarks for individual colors")

    # Misc
    parser.add_argument("-s", "--save", type=str, default=None,
                        help="filename to save figure")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="optional path to output processed seismograms")
    parser.add_argument("--noshow", action="store_true", default=False,
                        help="dont show the figure, default behavior will show")

    return parser.parse_args()


def find_nearest(array, value):
    """
    Find the nearest index in a consecutive (time) array given a chosen value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def _nearest_pow_2(x):
    """
    Find power of two nearest to x

    >>> _nearest_pow_2(3)
    2.0
    >>> _nearest_pow_2(15)
    16.0

    :type x: float
    :param x: Number
    :rtype: int
    :return: Nearest power of 2 to x
    """
    a = math.pow(2, math.ceil(np.log2(x)))
    b = math.pow(2, math.floor(np.log2(x)))
    if abs(a - x) < abs(b - x):
        return a
    else:
        return b


def convert_timezone(code, st):
    """
    When plotting in absolute time (args.time == "a*"), we allow time shifting 
    by time zone to get to the correct time local time. Returns streams with 
    converted time

    .. note::

        Confusingly, trying to input time zones in UTCDateTime like
        UTCDateTime("YYYY-MM-DDTHH:MM:SS-09) 
        assumes that the given date is in local time, and that the modifier -09
        is the conversion to UTC. But since seismic data is IN UTC, we cannot
        use this because it shifts us in the wrong way.

    :type code: str
    :param code: e.g., +08 to shift forward by 8 hours
    """
    assert(len(code) == 3), f"must be +?? or -??, not {code}"
    assert(code[0] in ["+", "-"])
    starttime = st[0].stats.starttime  # Dropping the 'Z' repr. UTC

    sign = code[0]
    shift = int(code[1:]) * 60 * 60  # hours

    if sign == "+":
        starttime += shift
    elif sign == "-":
        starttime -= shift

    for tr in st:
        if sign == "+":
            tr.stats.starttime += shift
        elif sign == "-":
            tr.stats.starttime -= shift

    return st


def spectrogram(ax, xvals, tr, per_lap=0.9, wlen=None, log=False,
                dbscale=False, mult=8., cmap=None, ncolors=256, zorder=None, 
                clip=[0.0, 1.0]):
    """
    Computes and plots spectrogram of the input data.

    .. note::
        Modified from ObsPys Spectrogram function
        https://docs.obspy.org/_modules/obspy/imaging/spectrogram.html\
            #spectrogram
    
        I also removed the non-log plotting option which uses imshow because I 
        didn't want to figure out how to change the x-axis values to match
        the waveform for sharex=True to work. So this might be a little slower
        when using non-log spectrogram y-scale but hopefully not significantly.


    :type tr: obspy Trace
    :param tr: Trace object to get data and sampling rate from
    :type per_lap: float
    :param per_lap: Percentage of overlap of sliding window, ranging from 0
        to 1. High overlaps take a long time to compute.
    :type wlen: int or float
    :param wlen: Window length for fft in seconds. If this parameter is too
        small, the calculation will take forever. If None, it defaults to a
        window length matching 128 samples.
    :type log: bool
    :param log: Logarithmic frequency axis if True, linear frequency axis
        otherwise.
    :type outfile: str
    :param outfile: String for the filename of output file, if None
        interactive plotting is activated.
    :type fmt: str
    :param fmt: Format of image to save
    :type axes: :class:`matplotlib.axes.Axes`
    :param axes: Plot into given axes, this deactivates the fmt and
        outfile option.
    :type dbscale: bool
    :param dbscale: If True 10 * log10 of color values is taken, if False the
        sqrt is taken.
    :type mult: float
    :param mult: Pad zeros to length mult * wlen. This will make the
        spectrogram smoother.
    :type cmap: :class:`matplotlib.colors.Colormap`
    :param cmap: Specify a custom colormap instance. If not specified, then the
        default ObsPy sequential colormap is used.
    :type zorder: float
    :param zorder: Specify the zorder of the plot. Only of importance if other
        plots in the same axes are executed.
    :type title: str
    :param title: Set the plot title
    :type show: bool
    :param show: Do not call `plt.show()` at end of routine. That way, further
        modifications can be done to the figure before showing it.
    :type clip: [float, float]
    :param clip: adjust colormap to clip at lower and/or upper end. The given
        percentages of the amplitude range (linear or logarithmic depending
        on option `dbscale`) are clipped.
    """
    if not cmap:
        cmap = obspy_sequential
    else:
        cmap = plt.get_cmap(cmap, ncolors)

    # enforce float for samp_rate
    samp_rate = float(tr.stats.sampling_rate)

    # set wlen from samp_rate if not specified otherwise
    if not wlen:
        wlen = 128 / samp_rate

    npts = tr.stats.npts

    # nfft needs to be an integer, otherwise a deprecation will be raised
    # XXX add condition for too many windows => calculation takes for ever
    nfft = int(_nearest_pow_2(wlen * samp_rate))

    if npts < nfft:
        msg = (f'Input signal too short ({npts} samples, window length '
               f'{wlen} seconds, nfft {nfft} samples, sampling rate '
               f'{samp_rate} Hz)')
        raise ValueError(msg)

    if mult is not None:
        mult = int(_nearest_pow_2(mult))
        mult = mult * nfft
    nlap = int(nfft * float(per_lap))

    data = tr.data - tr.data.mean()
    end = npts / samp_rate

    # Here we call not plt.specgram as this already produces a plot
    # specgram should be faster as it computes only the
    # arrays
    specgram, freq, time = mlab.specgram(data, Fs=samp_rate, NFFT=nfft,
                                         pad_to=mult, noverlap=nlap)

    if len(time) < 2:
        msg = (f'Input signal too short ({npts} samples, window length '
               f'{wlen} seconds, nfft {nfft} samples, {nlap} samples window '
               f'overlap, sampling rate {samp_rate} Hz)')
        raise ValueError(msg)

    # db scale and remove zero/offset for amplitude
    if dbscale:
        specgram = 10 * np.log10(specgram[1:, :])
    else:
        specgram = np.sqrt(specgram[1:, :])
    freq = freq[1:]

    # Normalize spectrogram so that lowest signal value is 0
    # !!! BC Is this okay?
    specgram -= specgram.min()

    vmin, vmax = clip
    if vmin < 0 or vmax > 1 or vmin >= vmax:
        msg = "Invalid parameters for clip option."
        raise ValueError(msg)
    _range = float(specgram.max() - specgram.min())
    vmin = specgram.min() + vmin * _range
    vmax = specgram.min() + vmax * _range
    norm = Normalize(vmin, vmax, clip=True)

    # calculate half bin width
    halfbin_time = (time[1] - time[0]) / 2.0
    halfbin_freq = (freq[1] - freq[0]) / 2.0

    kwargs = {'cmap': cmap, 'zorder': zorder}
    # pcolor expects one bin more at the right end, and we must also shift
    # our frequency bins to the center, currently they define the right edge
    freq = np.concatenate((freq, [freq[-1] + 2 * halfbin_freq]))
    freq -= halfbin_freq

    # We do the same for the time axis (NOTE: This is not used, see below)
    time = np.concatenate((time, [time[-1] + 2 * halfbin_time]))
    time -= halfbin_time

    # Hijack and change the time axis to match the waveforms so we can use
    # sharex, otherwise they are on separate x axes. 
    # We apply the same time shifting as above to get proper t0
    new_time = np.linspace(xvals.min(), xvals.max(), len(time))
    new_time += (time[0] + 2 * halfbin_time) / SECONDS_PER_DAY 

    # Plot times
    im = ax.pcolormesh(new_time, freq, specgram, norm=norm, **kwargs)

    # Log scaling for frequency values (y-axis)
    if log:
        ax.set_yscale('log')

    return im


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


def _set_xaxis_obspy_dates(ax, ticklabels_small=True, minticks=3, maxticks=6):
    """
    Set Formatter/Locator of x-Axis to use ObsPyAutoDateFormatter and do some
    other tweaking.

    In contrast to normal matplotlib ``AutoDateFormatter`` e.g. shows full
    timestamp on first tick when zoomed in so far that matplotlib would only
    show hours or minutes on all ticks (making it impossible to tell the date
    from the axis labels) and also shows full timestamp in matplotlib figures
    info line (mouse-over info of current cursor position).

    :type ax: :class:`matplotlib.axes.Axes`
    :rtype: None
    """
    ax.xaxis_date()
    locator = AutoDateLocator(minticks=minticks, maxticks=maxticks)
    locator.intervald[MINUTELY] = [1, 2, 5, 10, 15, 30]
    locator.intervald[SECONDLY] = [1, 2, 5, 10, 15, 30]
    ax.xaxis.set_major_formatter(ObsPyAutoDateFormatter(locator))
    ax.xaxis.set_major_locator(locator)
    if ticklabels_small:
        plt.setp(ax.get_xticklabels(), fontsize='small')


if __name__ == "__main__":
    args = parse_args()

    if not args.fid:
        sys.exit("positional argument `fid` required")

    # Populate Stream object
    st = Stream()
    for fid in args.fid:
        try:
            st += read(fid)
        except TypeError:
            st += read_sem(fid)
        print(fid)
    st.merge()

    # ==========================================================================
    #                           PROCESS WAVEFORMS
    # ==========================================================================
    if args.resample:
        print(f"resampling to {args.resample}")
        st.resample(sampling_rate=args.resample)

    # Time shift by requested amount
    if args.time.startswith("a") and len(args.time) > 1:
        st = convert_timezone(code=args.time[1:], st=st)

    # Trim data shorter so we don't process the entire waveform but add some
    # buffer so that preprocessing on the tails of the data doesn't show up.
    # the tail data will stay there but we will set the xlim of the plot to not
    # show it
    if args.xlim:    
        if args.time.startswith("a"):
            start = UTCDateTime(args.xlim[0])
            end = UTCDateTime(args.xlim[1])
            buffer = (end - start) * args.trim_pct  # some percentage of the record
            st.trim(starttime=start - buffer, endtime=end + buffer)
        else:
            starttime = st[0].stats.starttime
            buffer = (args.xlim[1] - args.xlim[0]) * args.trim_pct
            st.trim(starttime=starttime + float(args.xlim[0]) - buffer,
                    endtime=starttime + float(args.xlim[1]) + buffer)
        print(f"trimming with {buffer}s buffer on either end")

        if not st:
            print("Trimming removed all data, please check `xlim` values")
            sys.exit()

    # Preprocess waveforms
    taper = args.taper
    if args.integrate or args.differentiate:
        if args.taper == 0:
            taper = 0.05
            print(f"setting taper pct {taper}")
    if taper:
        st.taper(args.taper)
        print(f"tapering trace {taper * 100}%")
    if args.integrate:
        for i in range(args.integrate):
            print("integrating trace")
            st.integrate()
    if args.differentiate:
        for i in range(args.differentiate):
            print("differentiating trace")
            st.differentiate()

    # Allow different filters depending on min and max values given
    if args.fmin and args.fmax:
        print(f"bandpass {args.fmin}-{args.fmax}hz")
        st.filter("bandpass", freqmin=args.fmin, freqmax=args.fmax, 
                  zerophase=args.zerophase, corners=args.corners)
    elif args.fmin and not args.fmax:
        print(f"highpass {args.fmin}hz")
        st.filter("highpass", freq=args.fmin, zerophase=args.zerophase,
                  corners=args.corners)
    elif args.fmax and args.fmin is None:
        print(f"lowpass {args.fmax}hz")
        st.filter("lowpass", freq=args.fmax, zerophase=args.zerophase,
                  corners=args.corners)
    if args.fmin or args.fmax:
        print(f"zerophase={args.zerophase}")
        print(f"corners={args.corners}")

    # ==========================================================================
    #                           SET UP FIGURE
    # ==========================================================================
    if not args.spectrogram:
        f, ax = plt.subplots(figsize=(8, 4), dpi=200)
        axs = [ax]  # To play nice with some loops
        ax_spectra = None
    else:
        f, axs = plt.subplots(2, dpi=200, figsize=(8, 6), sharex=True)
        f.subplots_adjust(hspace=0)
        ax_spectra, ax = axs  # waveform on the bottom

    # ==========================================================================
    #                           PLOT WAVEFORM
    # ==========================================================================
    # Custom time axis
    if args.time.startswith("a"):
        # Set xvalues to datetime objects
        xvals = ((st[0].times() / SECONDS_PER_DAY) +
                        date2num(st[0].stats.starttime.datetime))
        _set_xaxis_obspy_dates(ax, maxticks=args.maxticks)
    else:
        xvals = st[0].times() 

        # Set time axis
        if args.time == "s":
            xvals /= 1  # not necessary but for consistency
        elif args.time == "m":
            xvals /= 60
        elif args.time == "h": 
            xvals /= 60 ** 2
        else:
            print("unknown time axis choice, default to 's'econds")
            xvals /= 1

        # Offset time axis based on user defined criteria
        xvals -= args.t0
        xvals += args.tstart
   
    for i, tr in enumerate(st):
        # Input a list of colors
        if len(args.colors) > 1:
            c = args.colors[i]
        # Input only a single color
        elif len(args.colors) == 1:
            # Allow C coloring
            if args.colors[0] == "C?":
                c = f"C{i}"
            else:
                c = args.colors[0]
        if args.labels:
            l = args.labels[i]
        else:
            l = tr.get_id()
        ax.plot(xvals, tr.data, c=c, lw=args.linewidth, zorder=6+i, label=l)

    # ==========================================================================
    #                  PLOT STREAM GAUGE (WARNING: SUPER CUSTOM)
    # ==========================================================================
    if args.stream_gage:
        assert args.time == "a-08", f"currently only works in AK local"

        # Read data from text file
        path = ("/Users/chow/Work/research/gulkanaseis24/data/USGS_data/"
                "phelan_creek_stream_guage_2024-09-07_to_2024-09-14.txt")
        assert(os.path.exists(path))

        data = np.loadtxt(path, skiprows=28, usecols=[2,4], delimiter="\t", 
                          dtype=str)
        times, height_ft = data.T  # time in AK local

        # Time is already in AK Local so we don't need to shift. If we did have
        # to then we would need to convert to UTC then shift by user request
        times = np.array([date2num(UTCDateTime(_).datetime) for _ in times])
        height_ft = np.array(height_ft, dtype=float)

        # Convert units
        if args.sg_units == "ft":
            height = height_ft
        elif args.sg_units == "in":
            height = height_ft * 12
        elif args.sg_units == "m":
            height = np.array([_ * 0.3048 for _ in height_ft.astype(float)])
        elif args.sg_units == "cm":
            height = np.array([_ * 30.48 for _ in height_ft.astype(float)])
        else:
            print("stream gage units `sg_units` should be in: ft, in, m, cm")
            sys.exit()
        
        # Subset data where we are plotting waveforms to get the correct ylims
        idx = np.where((times > xvals.min()) & (times < xvals.max()))

        if args.sg_rel:
            height -= height[idx].min()

        # Plot on the same axis as the waveform
        twax = ax.twinx()

        # # Plot the gradient, not very helpful
        # dx = (times[1] - times[0]) * 86400
        # grad_height = np.gradient(height[idx], dx)
        # twax.plot(times[idx], grad_height, "o-", lw=1, c="C0", 
        #         label="Phelan Cr. Gage", zorder=5, markersize=1.25, 
        #         alpha=0.5)

        twax.plot(times[idx], height[idx], "o-", lw=1, c="C0", 
                  label="Phelan Cr. Gage", zorder=5, markersize=1.25, 
                  alpha=0.5)

        _ylabel = f"Stream Height [{args.sg_units}]"
        if args.sg_rel:
            _ylabel = f"Relative {_ylabel}"
        twax.set_ylabel(_ylabel, rotation=-90, labelpad=20)

    # ==========================================================================
    #                           PLOT SPECTROGRAM
    # ==========================================================================
    # Don't plot spectrogram if we're plotting full data because it's too much
    if args.spectrogram and args.xlim:
        im = spectrogram(ax=ax_spectra, xvals=xvals, tr=st[args.sp_idx], 
                         log=args.sp_logscale, dbscale=args.sp_dbscale, 
                         cmap=args.sp_cmap, ncolors=args.sp_numcol
                         ) 
        ax_spectra.set_ylabel("Freq. [Hz]")
        ax_spectra.axis("tight")
        ax_spectra.grid(False)

        # SPECTROGRAM COLORBAR
        #  Pushes over the spectrogram in the same axis 
        # to make way for the colorbar
        _size = 2.5
        _pad = 0.05
        div_spectra = make_axes_locatable(ax_spectra)
        cax_spectra = div_spectra.append_axes("right", size=f"{_size}%", 
                                              pad=_pad)
        cbar = f.colorbar(im, cax=cax_spectra, orientation="vertical") 

        # Change the plot aeshtetic of the colorbar, should actually be in
        # set_plot_aesthetic() but I'm lazy
        for spine in cbar.ax.spines.values():
            spine.set_linewidth(1.5)
        cbar.ax.set_frame_on(True)
        # Padding was determined by trial and error
        if args.sp_dbscale:
            # _label = r"Rel. PSD [$10\log_{10}((m/s^2)/Hz)$]"
            _label = f"Rel. Power [dB]"

            _labelpad = -27
        else:
            _label = r"Rel. PSD [$(m/s^2/Hz)$]"
            _labelpad = -25
        cbar.ax.set_ylabel(_label, rotation=270, labelpad=_labelpad, fontsize=8)

        # HACKY: pushes over the waveform plot over by the same amount but '
        # then turns the space invisible to preserve the shared x-axis
        div = make_axes_locatable(ax)
        cax = div.append_axes("right", size=f"{_size}%", pad=_pad)
        cax.set_axis_off()

        # Have to do it for twax, too
        if args.stream_gage:
            div = make_axes_locatable(twax)
            cax = div.append_axes("right", size=f"{_size}%", pad=_pad)
            cax.set_axis_off()

    # ==========================================================================
    #                           PLOT TAUP ARRIVALS
    # ==========================================================================
    # Get phase arrivals from TauP if requested
    arrivals = None
    if args.tp_phases:
        assert(args.tp_dist is not None)
        assert(args.tp_depth is not None)
        dist_deg = kilometers2degrees(args.tp_dist)
        model = TauPyModel(model=args.tp_model)
        tp_arrivals = model.get_travel_times(source_depth_in_km=args.tp_depth,
                                             distance_in_degree=dist_deg,
                                             phase_list=args.tp_phases)

        # If some arrivals have multiple entires, only take first and last to 
        # get a range which we will plot with a window
        arrivals = {arrival.name: [] for arrival in tp_arrivals}
        for i, arrival in enumerate(tp_arrivals):
            arrivals[arrival.name].append(arrival.time)
            print(f"{arrival.name} = {arrival.time} s")

        if not arrivals:
            print(f"No arrivals found for given depth={args.tp_depth}km and "
                  f"distance {dist_deg:.2f}deg")
            
        arrival_dict = {}        
        for i, (name, times) in enumerate(arrivals.items()):
            if times[0] == times[-1]:
                alpha = 1
            else: 
                alpha = 0.3

            # Figure out the maximum amplitude in this time window 
            win_start = find_nearest(xvals, times[0])
            win_end = find_nearest(xvals, times[-1]) + 1
            max_amp = np.amax(st[0].data[win_start:win_end])
            arrival_dict[name] = max_amp

            plt.axvspan(times[0], times[-1], label=f"{name} ({max_amp:.2E})", 
                        color=f"C{i}", alpha=alpha, zorder=7)
        plt.legend(fontsize=8, loc="upper left", frameon=False)

    # ==========================================================================
    #                           PLOT TMARKS
    # ==========================================================================
    if args.tmarks:
        if len(args.tmarks_c) == 1:
            colors = args.tmarks_c * len(args.tmarks)
        else:
            colors = args.tmarks_c
        for tmark, c in zip(args.tmarks, colors):
            if args.time.startswith("a"):
                tmark = date2num(UTCDateTime(tmark).datetime)
                # date2num(st[0].stats.starttime.datetime) + tmark
            ax.axvline(tmark, c=c, lw=0.5)

    # ==========================================================================
    #                           PLOT AESTHETICS
    # ==========================================================================
    f.legend(loc="upper right", bbox_to_anchor=(1,1), 
             bbox_transform=ax.transAxes, prop={"size": 6})

    if args.time.startswith("a"):
        ax.set_xlabel(f"Time [UTC{args.time[1:]}]")
    else:
        ax.set_xlabel(f"Time [{args.time}]")
    ax.set_ylabel(args.ylabel)

    # Subset x axis
    xstart, xend = xvals.min(), xvals.max()
    if args.xlim:
        if args.time.startswith("a"):
            xstart = date2num(UTCDateTime(args.xlim[0]).datetime)
            xend = date2num(UTCDateTime(args.xlim[1]).datetime)
        else:
            xstart, xend = [float(_) for _ in args.xlim]
    ax.set_xlim(xstart, xend)

    # Subset y axis for waveform plot
    ymax = np.amax([st[0].data.min(), st[0].data.max()])
    ymin = -1 * ymax
    if args.ylim:
        # Allow for one entry to set min/max if they're the same
        if len(args.ylim) == 1:
            ymin, ymax = [-1 * args.ylim[0], args.ylim[0]]
        else:
            ymin, ymax = args.ylim
    else:
        ymax = np.amax([st[0].data.min(), st[0].data.max()])
        ymin = -1 * ymax
    
    ax.set_ylim(ymin, ymax)

    # Spectrogram annotation if there are multiple waveforms plotted
    if args.spectrogram and len(st) > 1:
        ax_spectra.text(0.99, 0.99, st[args.sp_idx].get_id(), 
                        horizontalalignment="right", verticalalignment="top", 
                        transform=ax_spectra.transAxes, fontsize=8)

    # Finish off by setting plot aesthetics
    if args.title is None:
        title = f"{st[0].stats.starttime.year}.{st[0].stats.starttime.julday}"
        if st[0].stats.starttime.julday != st[0].stats.endtime.julday:
            title += f"-{st[0].stats.endtime.julday}"

        if args.fmin or args.fmax:
            _fmin = args.fmin or 0
            _fmax = args.fmax or st[0].stats.sampling_rate / 2
            title += f" [{_fmin}, {_fmax}]Hz"

        # Append some information on the TauP arrivals
        if arrivals:
            title += (f"\n(TauP={args.tp_model}; $\\Delta$={args.tp_dist}km; "
                      f"Z={args.tp_depth}km)")
        if args.title_append:
            title += f"\n{args.title_append}"
    else:
        title = args.title
    plt.suptitle(title)

    # Final plotting touches
    if args.labels:
        plt.legend()
    set_plot_aesthetic(ax)
    if ax_spectra:
        set_plot_aesthetic(ax_spectra, ytick_format="plain")
    f.tight_layout()

    # ==========================================================================
    #                           FINALIZE PLOT
    # ==========================================================================
    _transparent = False  # toggle for transparency in the background
    if args.save == "auto":
        plt.savefig(f"{args.fid[0]}.png", transparent=_transparent)
    elif args.save is not None:
        plt.savefig(args.save, transparent=_transparent)

    if not args.noshow:
        plt.show()
    plt.close("all")

    # ==========================================================================
    #                               OUTPUT
    # ==========================================================================
    if args.output:
        if not os.path.exists(args.output):
            os.makedirs(args.output)
        st.write(os.path.join(args.output, args.fid), format="MSEED")

