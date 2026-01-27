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
import matplotlib as mpl
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

# Only needed for SPECFEM synthetics, ignore if not installed
try:
    from pysep import read_sem
except ImportError:
    read_sem = None
    pass

SECONDS_PER_DAY = 3600.0 * 24.0

def parse_args():
    """All modifications are accomplished with command line arguments"""
    parser = argparse.ArgumentParser()

    # Waveform Processing
    parser.add_argument("fids", nargs="+", 
                        help="required, file ID(s), wildcards accepted")
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
    parser.add_argument("--detrend", nargs="?", type=str, default=None,
                        help="demean the data prior to other processing",
                        choices=["simple", "linear", "demean"])
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
    
    # Waveform options
    parser.add_argument("--wf_abs", action="store_true", 
                        help="plot the absolute value of the waveform, does" 
                             "not affect processing, just final plotting")
    parser.add_argument("--wf_type", nargs="?", type=str, default="default",
                        help="Option for how to plot the waveforms:\n"
                             "'default': plot all `fids` on top of each other " 
                             "with absolute amplitudes.\n"
                             "'stack': stack all `fids` on top of each other w " 
                             "some transparency, and then plot the mean value " 
                             "of all traces as a dark line on top of all\n"
                             "'recsec': plot all waveforms with relative " 
                             "amplitudes on a record section style plot so that" 
                             "each waveform is visually distinct from the " 
                             "others\n"
                             "'recsec_stack', recsec but add stack",
                        choices=["default", "stack", "recsec", "recsec_stack"]
                             )
    parser.add_argument("--wf_norm", action="store_true", 
                        help="normalize traces to their individual max")
    parser.add_argument("--wf_recsec_spacing", nargs="?", type=float, default=1,
                        help="if `wf_type` is some form of recsec, choose the " 
                             "spacing modifier between each of the waveforms. " 
                             "By default it matches the order of magnitude so " 
                             "this modifier set to 1 is a good starting guess "
                             "but you can move it up or down by steps of 0.1 " 
                             "to get waveforms closer or further")
    parser.add_argument("--wf_spacing_exact", action="store_true", 
                        help="For `wf_recsec_spacing` the value is used exact "
                             "rather than as an order of magnitude scale")
    parser.add_argument("--wf_order", type=int, nargs="+", default=None,
                        help="allow custom order when plotting waveforms, " \
                             "based on the alphabetical index. First you need "
                             "to plot the normal sorted one and then you can" \
                             "rearrange based on that order")
    
    # Plot Aesthetics for Waveforms
    parser.add_argument("-c", "--colors", nargs="+", type=str, default="k",
                        help="how to select color for waveforms:\n" 
                             "'k': select one color for all waveforms\n" 
                             "['k', 'r']: make a list of colors for each trace "
                             "must match the length of traces\n" 
                             "'CN': use CN colors, predefiend colors in MPL,\n" 
                             "'viridis': use a colormap and plot discrete " 
                             "colors from the map based on length of stream"
                             )
    parser.add_argument("--alphas", nargs="+", type=float, default=None,
                        help="alpha of the time series line, number of inputs "
                             "must match the length of `fid`")
    parser.add_argument("-l", "--labels", nargs="+", type=str, default=None,
                        help="optional labels legend, must match len of `fid`")
    parser.add_argument("-lw", "--linewidth", nargs="?", type=float, 
                        default=0.5, help="linewidth of the time series line")
    parser.add_argument("--ylabel", nargs="?", type=str, default="Amplitude",
                        help="label for units, defaults to displacement")
    parser.add_argument("-y", "--ylim", nargs="+", type=float, default=None,
                        help="amplitude axis limits in s")
    
    # Time axis (X-axis)
    parser.add_argument("-t", "--time", nargs="?", type=str, default="s",
                    help="units for x-axis/time axis. choice: 's'econds "
                            "(default), 'm'inutes, 'h'ours, 'a'bsolute (wip)."
                            "If using 'a' you may add '+i' or '-i' to "
                            "time shift the array, e.g., to go from UTC to "
                            "local time. E.g., 'a-7' will subtract 7 hours.")
    parser.add_argument("--minticks", type=int, default=3, 
                        help="min ticks if --time='a'")
    parser.add_argument("--maxticks", type=int, default=6, 
                        help="max ticks if --time='a'")
    parser.add_argument("-x", "--xlim", nargs="+", default=None,
                        help="time axis limits in s or if `time`=='a' then "
                             "values should be in datetime, see tmarks")
    parser.add_argument("-tm", "--tmarks", nargs="+", 
                        help="plot vertical lines at given relative times, "
                             "should match the units of `time`. If `time`=='a' "
                             "then each tmark should be a datetime "
                             "YYYY-MM-DDTHH:MM:SS")
    parser.add_argument("--tmarks_c", nargs="+", default="k",
                        help="colors for each of the time marks, should "
                             "either be single letter for all marks or match " 
                             "length of tmarks for individual colors")

    # TauP Phase Arrivals
    parser.add_argument("--tp_phases", nargs="+", type=str, default=None,
                        help="taup phase arrivals, requires all 'tp_*' params")
    parser.add_argument("--tp_model", nargs="?", type=str, default="iasp91",
                        help="taup model, defaults to 'iasp91'")
    parser.add_argument("--tp_dist_km", nargs="?", type=float, default=None,
                        help="TauP source receiver distance in km")
    parser.add_argument("--tp_dist_deg", nargs="?", type=float, default=None,
                        help="TauP source receiver distance in degrees")
    parser.add_argument("--tp_depth", nargs="?", type=float, default=None,
                        help="TauP source depth km")
    parser.add_argument("--tp_start", nargs="?", type=str, default=None,
                        help="origintime of the event, must match the format "
                             "of --time. If not given, assumes first sample "
                             "of --xlim defines origintime")
    
    # Spectrogram (all parameters end with _s)
    parser.add_argument("--spectrogram", action="store_true", default=False,
                        help="plot spectrogram of the raw trace. See all " 
                             "'sp_*' parameters to control the look")
    parser.add_argument("--sp_cmap", nargs="?", type=str, 
                        default="viridis",
                        help="colormap of the spectrogram")
    parser.add_argument("--sp_numcol", nargs="?", type=int, default=256,
                        help="number of colors in colormap of the spectrogram")
    parser.add_argument("--sp_dbscale_off", action="store_true",
                        help="Turn off dB scaling for colors/amplitudes in " 
                             "spectrogram")    
    parser.add_argument("--sp_logscale", action="store_true",
                        help="turn on log scale for spectrogram y-axis")
    parser.add_argument("--sp_clip", default=[0., 1.], type=float, nargs="+",
                        help=" adjust cmap to clip at lower and/or upper end")
    parser.add_argument("--sp_vmax", default=None, type=float,
                        help="maximum value threshold for spectrogram")
    parser.add_argument("--sp_idx", type=int, default=0,
                        help="Iff multiple waveforms plotted, choose which of "
                             "them to use for the spectrogram. Defaults to 0. " 
                             "order is alphabetical")
    
    # Plot additional time series
    parser.add_argument("--add_trace", nargs="?", type=str, default=None,
                        help="Path to plot an additional time series that does "
                             "not have the same sampling rate . "
                             "All --tr_* pars are associated. Data should be " 
                             "in 2-column ASCII or an ObsPy readable format")
    parser.add_argument("--tr_time", nargs="?", type=str, default="a",
                        help="Time of the time series [a]bsolute (UTC) or " \
                             "[r]elative to xlim[0] or start of main trace",
                        choices=["a", "r"])
    parser.add_argument("--tr_ylabel", type=str, default="",
                        help="Y-label text")
    parser.add_argument("--tr_label", type=str, default="",
                        help="Label text for shared legend")
    parser.add_argument("--tr_color", type=str, default="C0",
                       help="color of the time series")
    
    # CUSTOM
    parser.add_argument("--stream_gage", action="store_true", default=False,
                        help="For GULKANASEIS data only, plots stream gauge " 
                             "data at the bottom of the waveform plot with a " 
                             "twin X axis")

    # Misc plotting options
    parser.add_argument("--fig_size", type=float, nargs="+", default=None,
                        help="Figure size")
    parser.add_argument("--dpi", type=float, nargs="?", default=None,
                        help="Dots per inch, default 200")
    parser.add_argument("--fig_len", type=float, nargs="?", default=8,
                        help="Figure length, use with 'fig_asp', overrides "
                             "'--fig_size' if both given")
    parser.add_argument("--fig_asp", type=float, nargs="?", default=None,
                        help="Figure aspect ratio, use with 'fig_len', "
                             "overrides '--fig_size' if given")
    parser.add_argument("--no_legend", action="store_true",
                        help="Turn off waveform legend")
    parser.add_argument("--ncol_legend", type=int, default=1,
                        help="Number of columns in legend")
    parser.add_argument("--title", nargs="?", type=str, default=None,
                        help="title of the figure, defaults to ID and fmin/max")
    parser.add_argument("-ta", "--title_append", nargs="?", type=str, 
                        default="", help="append to default title")
    parser.add_argument("--frameoff", action="store_true",
                        help="Turn everything off except the waveform")
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


def convert_timezone(code, st):
    """
    When plotting in absolute time (args.time == "a*"), we allow time 
    shifting  by time zone to get to the correct time local time. Returns 
    streams with converted time

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


class PrettyPlot():
    """
    Command line and scriptable waveform and spectrogram plotter
    """
    def __init__(self, fids, 
                 # Processing
                 taper=0, fmin=None, fmax=None, zerophase=False, corners=4,
                 resample=False, t0=0, tstart=0, detrend=False, integrate=0,
                 differentiate=0, trim_pct=1.,
                 # Waveform plotting
                 wf_abs=False, wf_type="default", wf_norm=False, 
                 wf_recsec_spacing=1, wf_spacing_exact=False, wf_order=None,
                 # Plotting Aesthetics
                 colors="k", alphas=None, labels=None, linewidth=0.5, 
                 ylabel="amplitude", ylim=None,
                 # Time Axis
                 time="s", minticks=3, maxticks=6, xlim=6, tmarks=None, 
                 tmarks_c="k",
                 # TauP 
                 tp_phases=None, tp_model="iasp91", tp_dist_km=None,
                 tp_dist_deg=None, tp_depth=None, tp_start=None,
                 # Spectrogram
                 spectrogram=False, sp_cmap="viridis", sp_numcol=256,
                 sp_dbscale=True, sp_logscale=False, sp_clip=None,
                 sp_vmax=None, sp_idx=0,
                 # Additional Time Series
                 add_trace=None, tr_time="a", tr_label="", tr_ylabel="", 
                 tr_color="C0", 
                 stream_gage=False,
                 # Misc.
                 fig_size=None, dpi=None, fig_len=None, fig_asp=None, 
                 legend=True, ncol_legend=1, title=None, title_append="",
                 save=None, output=None, show=True, frameoff=False,
                 **kwargs
                 ):
        """Input parameters, see argparser for descriptions"""
        self.fids = fids
        self.taper = taper
        self.fmin = fmin
        self.fmax = fmax
        self.zerophase = zerophase
        self.corners = corners
        self.resample = resample
        self.t0 = t0
        self.tstart = tstart
        self.detrend = detrend
        self.integrate = integrate
        self.differentiate = differentiate
        self.trim_pct = trim_pct

        self.wf_abs = wf_abs
        self.wf_type = wf_type
        self.wf_norm = wf_norm
        self.wf_recsec_spacing = wf_recsec_spacing
        self.wf_spacing_exact = wf_spacing_exact
        self.wf_order = wf_order

        self.colors = colors
        self.alphas = alphas
        self.labels = labels
        self.linewidth = linewidth
        self.ylim = ylim
        self.ylabel = ylabel

        self.time = time
        self.maxticks = maxticks
        self.minticks = minticks
        self.xlim = xlim
        self.tmarks = tmarks
        self.tmarks_c = tmarks_c

        self.tp_phases = tp_phases
        self.tp_model = tp_model
        self.tp_dist_km = tp_dist_km
        self.tp_dist_deg = tp_dist_deg
        self.tp_depth = tp_depth
        self.tp_start = tp_start
        self.spectrogram = spectrogram
        self.sp_cmap = sp_cmap
        self.sp_numcol = sp_numcol
        self.sp_dbscale = sp_dbscale
        self.sp_logscale = sp_logscale
        self.sp_clip = sp_clip or [0., 1.]
        self.sp_vmax = sp_vmax
        self.sp_idx = sp_idx

        self.add_trace = add_trace
        self.tr_time = tr_time
        self.tr_label = tr_label
        self.tr_ylabel = tr_ylabel
        self.tr_color = tr_color

        self.fig_size = fig_size
        self.dpi = dpi
        self.fig_len = fig_len
        self.fig_asp = fig_asp
        self.legend = legend
        self.ncol_legend = ncol_legend
        self.title = title
        self.title_append = title_append
        self.save = save
        self.output = output
        self.show = show
        self.frameoff = frameoff

        # Populate Stream object
        self.st = Stream()
        for fid in self.fids:
            try:
                self.st += read(fid)
            except TypeError:
                if bool(read_sem):
                    self.st += read_sem(fid)
                else:
                    sys.exit("no function `read_sem()` from PySEP")
            print(fid)
        print(f"{self.st.__str__(extended=True)}")

    def setup_plot(self, dpi=200):
        """
        Set up the plot based on input parameters
        """
        if self.spectrogram:
            if self.fig_len and self.fig_asp:
                figsize = (self.fig_len, self.fig_len * self.fig_asp)
            else:
                figsize = self.fig_size or (8, 6)  # default

            dpi = self.dpi or dpi

            print("\tsetting up waveform and spectrogram plot")
            self.f, axs = plt.subplots(2, dpi=dpi, figsize=figsize, sharex=True)
            self.f.subplots_adjust(hspace=0)
            self.ax_spectra, self.ax = axs  # waveform on the bottom
        else:
            if self.fig_len and self.fig_asp:
                figsize = (self.fig_len, self.fig_len * self.fig_asp)
            else:
                figsize = self.fig_size  or (7, 4)  # default

            dpi = self.dpi or dpi

            print("\tsetting up waveform plot")
            self.f, self.ax = plt.subplots(dpi=dpi, figsize=figsize)

    def trim_waveform(self):
        """Trim the waveform data so we don't plot the whole thing"""         
        if self.resample:
            print(f"\tresampling to {self.resample}")
            self.st.resample(sampling_rate=self.resample)
       
        # Shift to correct time axis
        if self.time.startswith("a") and len(self.time) > 1:
            self.st = convert_timezone(code=self.time[1:], st=self.st)

        # Trim data shorter so we don't process the entire waveform but add some
        # buffer so that preprocessing on the tails of the data doesn't show up.
        # the tail data will stay there but we will set the xlim of the plot to 
        # not show it
        if self.xlim:    
            if self.time.startswith("a"):
                start = UTCDateTime(self.xlim[0])
                end = UTCDateTime(self.xlim[1])
                buffer = (end - start) * self.trim_pct  # some % of the record
                self.st.trim(starttime=start - buffer, endtime=end + buffer)
            else:
                starttime = self.st[0].stats.starttime
                buffer = (float(self.xlim[1]) - float(self.xlim[0])) * \
                                                                   self.trim_pct
                self.st.trim(
                    starttime=starttime + (float(self.xlim[0]) - buffer),
                    endtime=starttime + (float(self.xlim[1]) + buffer)
                    )
            print(f"\ttrimming with {buffer}s buffer on either end")

            if not self.st:
                print("Trimming removed all data, please check `xlim` values")
                sys.exit()

    def process_waveforms(self):
        """
        Process the waveforms based on input parameters
        """
        if self.detrend:
            print(f"\tdetrending with '{self.detrend}'")
            self.st.detrend(type=self.detrend)
        taper = self.taper
        if self.integrate or self.differentiate:
            if self.taper == 0:
                taper = 0.05
                print(f"\tsetting taper pct {taper}")
        if taper:
            self.st.taper(taper)
            print(f"\ttapering trace {taper * 100}%")
        if self.integrate:
            for i in range(self.integrate):
                print("integrating trace")
                self.st.integrate()
        if self.differentiate:
            for i in range(self.differentiate):
                print("differentiating trace")
                self.st.differentiate()
        if self.fmin and self.fmax:
            print(f"\tfiltering between {self.fmin} and {self.fmax} Hz, "
                  f"zerophase={self.zerophase}, corners={self.corners}")
            self.st.filter(
                "bandpass", freqmin=self.fmin, freqmax=self.fmax,
                zerophase=self.zerophase, corners=self.corners
                )
        elif self.fmin:
            print(f"\thighpass filtering at {self.fmin} Hz, "
                  f"zerophase={self.zerophase}, corners={self.corners}")
            self.st.filter(
                "highpass", freq=self.fmin, zerophase=self.zerophase,
                corners=self.corners
                )
        elif self.fmax:
            print(f"\tlowpass filtering at {self.fmax} Hz, "
                  f"zerophase={self.zerophase}, corners={self.corners}")
            self.st.filter(
                "lowpass", freq=self.fmax, zerophase=self.zerophase,
                corners=self.corners
                )
        
        self.st.merge()

    def plot_waveforms(self):
        """
        Plot the waveforms based on input parameters
        """
        if self.time.startswith("a"):
            # Set xvalues to datetime objects
            xvals = ((self.st[0].times() / SECONDS_PER_DAY) +
                            date2num(self.st[0].stats.starttime.datetime))
            _set_xaxis_obspy_dates(self.ax, minticks=self.minticks, 
                                   maxticks=self.maxticks)
        else:
            xvals = self.st[0].times() 

            # Set time axis
            if self.time == "s":
                xvals /= 1  # not necessary but for consistency
            elif self.time == "m":
                xvals /= 60
            elif self.time == "h": 
                xvals /= 60 ** 2
            else:
                print("unknown time axis choice, default to 's'econds")
                xvals /= 1

            # Offset time axis based on user defined criteria
            xvals -= self.t0
            xvals += self.tstart
        
        if not self.wf_order:
            data = np.array([tr.data for tr in self.st])
            labels = [tr.get_id() for tr in self.st]
        else:
            # Allow arbitrary order set by User
            data, labels = [], []
            for i in self.wf_order:
                data.append(self.st[i].data)
                labels.append(self.st[i].get_id())
        n = len(data)
        # Overwrite waveform labels for the legend if User requests
        if self.labels:
            labels = self.labels
        # Determine how many waveforms we will be plotting
        if self.wf_type.endswith("stack"):
            n += 1
        # Set up the type of waveform
        if self.wf_type == "stack":  # stack or recsec_stack
            print("waveform option `stack`")
            print("overriding `alphas` and `colors` for action `stack`")
            # Adjust these if you want 
            alphas = [0.5] * n
            colors = [f"C{i}" for i in range(n)]
        else:
            # Allow list of alpha values
            if not self.alphas:
                alphas = [1] * n
            else:
                alphas = self.alphas
            # Allow single color, list of colors, or C colors
            if len(self.colors) > 1:
                colors = self.colors
            elif len(self.colors) == 1:
                if self.colors[0] == "C?":
                    colors = [f"C{i}" for i in range(n)]
                # Allow colormaps that we make into discrete colors
                elif len(self.colors[0]) > 1:
                    cmap = mpl.colormaps[self.colors[0]]
                    colors = cmap(np.linspace(0.1, 1, n))
                else:
                    colors = [self.colors[0]] * n

        # Allow for absolute amplitudes
        if self.wf_abs:
            data = [np.abs(d) for d in data]
            labels = [f"abs {l}" for l in labels]
        # Generate the stacked data (by default it takes the mean)
        if self.wf_type.endswith("stack"):
            wf_stack_type = "absmean"
            # Determine how to stack
            if wf_stack_type == "mean":
                stacked_data = np.mean(data, axis=0)
            # Pseudo envelope by meaning all abs values
            elif wf_stack_type == "absmean":
                stacked_data = np.abs(np.mean(np.abs(data), axis=0))
            # Add stacked data to the list of plottable data
            data = np.vstack((data, stacked_data))
            labels.append("abs stacked mean")
        # Normalize the dat arrays
        if self.wf_norm:
            for i, d in enumerate(data[:]):
                data[i] = d / d.max()
        # Value shift the data array to make a record section
        if self.wf_type.startswith("recsec"):
            print("waveform option `recsec`, scaling y-axis")
            if self.wf_spacing_exact:
                # Directly set the spacing with a vlaue
                spacing = self.wf_recsec_spacing
            else:
                # Set an arbitrary spacing based on the order of magnitude 
                oom = np.floor(np.log10(np.amax(self.st[0].data)))
                spacing = self.wf_recsec_spacing* 10 ** oom
            print(f"\trecsec spacing is {spacing}")
            for i, d in enumerate(data[:]):
                data[i] = d + (i * spacing)
        # Plot the waveforms
        for i, d in enumerate(data):
            self.ax.plot(
                xvals, d, c=colors[i], lw=self.linewidth, zorder=6+i, 
                label=labels[i], alpha=alphas[i]
                )
            
        self._xvals = xvals

    def plot_stream_gage(self, relative=False, units="m"):
        """
        Experimental Phelan Creek Stream Gage
        """
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
        if units == "ft":
            height = height_ft
        elif units == "in":
            height = height_ft * 12
        elif units == "m":
            height = np.array([_ * 0.3048 for _ in height_ft.astype(float)])
        elif units == "cm":
            height = np.array([_ * 30.48 for _ in height_ft.astype(float)])
        else:
            print("stream gage units `sg_units` should be in: ft, in, m, cm")
            sys.exit()
        
        # Subset data where we are plotting waveforms to get the correct ylims
        idx = np.where((times > self.xvals.min()) & (times < self.xvals.max()))

        if relative:
            height -= height[idx].min()

        # Plot on the same axis as the waveform
        twax = self.ax.twinx()

        twax.plot(times[idx], height[idx], "o-", lw=1, c="C0", 
                label="Phelan Cr. Gage", zorder=5, markersize=1.25, 
                alpha=0.5)

        _ylabel = f"Stream Height [{units}]"
        if relative:
            _ylabel = f"Relative {_ylabel}"
        twax.set_ylabel(_ylabel, rotation=-90, labelpad=20)
            
    def plot_additional_traces(self):
        """Plot additional time series if requested """
        twax = self.ax.twinx()
        st_new = read(self.add_trace).merge()
        data_new = st_new[0].data

        # Custom time axis
        if self.tr_time.startswith("a"):
            # Set xvalues to datetime objects
            xvals_new = ((st_new[0].times() / SECONDS_PER_DAY) +
                            date2num(st_new[0].stats.starttime.datetime))
            _set_xaxis_obspy_dates(twax, maxticks=self.maxticks)
        else:
            xvals_new = st_new[0].times() 

            # Set time axis
            if self.tr_time == "s":
                xvals_new /= 1  # not necessary but for consistency
            elif self.tr_time == "m":
                xvals_new /= 60
            elif self.tr_time == "h": 
                xvals_new /= 60 ** 2
            else:
                print("unknown time axis choice, default to 's'econds")
                xvals_new /= 1

            # Offset time axis based on user defined criteria
            xvals_new -= self.t0
            xvals_new += self.tstart

        # Plot on the same axis as the waveform
        if self.tr_label:
            tr_label = self.tr_label
        else:
            tr_label = st_new[0].get_id()
        twax.plot(xvals_new, data_new, lw=1, c=self.tr_color, 
                  label=tr_label, zorder=5, markersize=1.25, alpha=0.4)

        twax.set_ylabel(self.tr_ylabel, rotation=-90, labelpad=20)

        if self.spectrogram:
            _size = 2.5
            _pad = 0.05
            div = make_axes_locatable(twax)
            cax = div.append_axes("right", size=f"{_size}%", pad=_pad)
            cax.set_axis_off()

    def plot_spectrogram(self):    
        """Plot spectrogram on a separate axis"""
        im = spectrogram(
            ax=self.ax_spectra, xvals=self._xvals, tr=self.st[self.sp_idx], 
            log=self.sp_logscale, dbscale=self.sp_dbscale, 
            cmap=self.sp_cmap, ncolors=self.sp_numcol,
            clip=self.sp_clip,
            ) 
        self.ax_spectra.set_ylabel("Freq. [Hz]")
        self.ax_spectra.axis("tight")
        self.ax_spectra.grid(False) 

        # SPECTROGRAM COLORBAR
        #  Pushes over the spectrogram in the same axis
        _size = 2.5
        _pad = 0.05
        div_spectra = make_axes_locatable(self.ax_spectra)
        cax_spectra = div_spectra.append_axes("right", size=f"{_size}%",
                                                pad=_pad)
        cbar = self.f.colorbar(im, cax=cax_spectra, orientation="vertical")
        # Change the plot aeshtetic of the colorbar, should actually be in
        # set_plot_aesthetic() but I'm lazy
        for spine in cbar.ax.spines.values():
            spine.set_linewidth(1.5)
        cbar.ax.set_frame_on(True)
        # Padding was determined by trial and error
        if self.sp_dbscale:
            # _label = r"Rel. PSD [$10\log_{10}((m/s^2)/Hz)$]"
            _label = f"Power [dB]"
            _labelpad = -37.5
        else:
            _label = r"PSD [$(m/s^2/Hz)$]"
            _labelpad = -25
        cbar.ax.set_ylabel(_label, rotation=270, labelpad=_labelpad, fontsize=8)

        # HACKY: pushes over the waveform plot over by the same amount but '
        # then turns the space invisible to preserve the shared x-axis
        div = make_axes_locatable(self.ax)
        cax = div.append_axes("right", size=f"{_size}%", pad=_pad)
        cax.set_axis_off()

    def plot_taup_arrivals(self):
        """Get and plot TauP arrivals if requested"""
        if not self.tp_phases:
            return

        # Get phase arrivals from TauP if requested
        arrivals = None
        assert(self.tp_dist_km is not None or self.tp_dist_deg is not None)
        assert(self.tp_depth is not None)
        if self.tp_dist_km:
            dist_deg = kilometers2degrees(self.tp_dist_km)
        else:
            dist_deg = self.tp_dist_deg

        model = TauPyModel(model=self.tp_model)
        tp_arrivals = model.get_travel_times(
            source_depth_in_km=self.tp_depth,
            distance_in_degree=dist_deg,
            phase_list=self.tp_phases
            )

        # If some arrivals have multiple entires, only take first and last to 
        # get a range which we will plot with a window
        arrivals = {arrival.name: [] for arrival in tp_arrivals}
        print(f"\tTauP Arrivals for {self.tp_model}")
        for i, arrival in enumerate(tp_arrivals):
            arrivals[arrival.name].append(arrival.time)
            print(f"\t\t{arrival.name} = {arrival.time:.2f} s")

        if not arrivals:
            print(f"\tNo arrivals found for given depth={self.tp_depth}km and "
                  f"distance {dist_deg:.2f}deg")
           
        # Determine the time series starttime of the event because the TauP
        # times are relative to an origin time
        if not self.tp_start: 
            if self.time.startswith("a"):
                tp_start = UTCDateTime(self.xlim[0])
            else:
                tp_start = float(self.xlim[0])
        else:
            if self.time.startswith("a"):
                tp_start = UTCDateTime(self.tp_start)
            else:
                tp_start = float(self.tp_start)  # seconds

        for i, (name, times) in enumerate(arrivals.items()):
            # Sometimes there are multiple arrival times for the same phase
            if times[0] == times[-1]:
                alpha = 1
            else: 
                alpha = 0.3

            # Convert arrival times to time series reference
            if self.time.startswith("a"):
                times = [date2num((tp_start + time).datetime) for time in times]
            else:
                times = [tp_start + time for time in times]
            # Figure out the maximum amplitude in this time window
            win_start = find_nearest(self._xvals, times[0])
            win_end = find_nearest(self._xvals, times[-1]) + 1
            max_amp = np.amax(self.st[0].data[win_start:win_end])
            self.ax.axvspan(
                times[0], times[-1], label=f"{name} ({max_amp:.2E})", 
                color=f"C{i}", alpha=alpha, zorder=7
                )
        # plt.legend(prop={"size": 2}, loc="upper left", frameon=False

    def plot_tmarks(self):
        """
        Plot time marks if requested
        """
        if len(self.tmarks_c) == 1:
            colors = self.tmarks_c * len(self.tmarks)
        else:
            colors = self.tmarks_c

        for tmark, c in zip(self.tmarks, colors):
            if self.time.startswith("a"):
                tmark = date2num(UTCDateTime(tmark).datetime)
                # date2num(st[0].stats.starttime.datetime) + tmark
            self.ax.axvline(tmark, c=c, lw=0.5)

    def plot_aesthetics(self):
        """
        Set plot aesthetics
        """
        if self.legend:
            self.ax.legend(loc="upper right", bbox_to_anchor=(1,1), 
                          bbox_transform=self.ax.transAxes, 
                          prop={"size": 5},
                          ncol=self.ncol_legend, fontsize='tiny')

        if self.time.startswith("a"):
            self.ax.set_xlabel(f"Time [UTC{self.time[1:]}]")
        else:
            self.ax.set_xlabel(f"Time [{self.time}]")
        if self.ylabel:
            self.ax.set_ylabel(self.ylabel)
        else:
            self.ax.set_ylabel("Amplitude")

        # Subset x axis
        xstart, xend = self._xvals.min(), self._xvals.max()
        if self.xlim:
            if self.time.startswith("a"):
                xstart = date2num(UTCDateTime(self.xlim[0]).datetime)
                xend = date2num(UTCDateTime(self.xlim[1]).datetime)
            else:
                xstart, xend = [float(_) for _ in self.xlim]
        self.ax.set_xlim(xstart, xend)  

        # Subset y axis for waveform plot
        if self.ylim:
            if len(self.ylim) == 2:
                ymin, ymax = self.ylim
            else:
                ymax = abs(self.ylim[0])
                ymin = -1 * ymax
            self.ax.set_ylim(ymin, ymax)

        # Spectrogram annotation if there are multiple waveforms plotted
        # put it outside the axis so it shows up on white background
        if self.spectrogram and len(self.st) > 1:
            self.ax_spectra.text(0.99, 1.01, self.st[self.sp_idx].get_id(), 
                            horizontalalignment="right", 
                            verticalalignment="bottom", 
                            transform=self.ax_spectra.transAxes, fontsize=8)

        # Finish off by setting plot aesthetics
        set_plot_aesthetic(self.ax)

        if self.spectrogram:
            set_plot_aesthetic(self.ax_spectra, ytick_format="plain")

        if self.title is None:
            title = f"{self.st[0].stats.starttime.year}."
            title += f"{self.st[0].stats.starttime.julday:0>3}"
            if self.st[0].stats.starttime.julday != \
                self.st[0].stats.endtime.julday:
                title += f"-{self.st[0].stats.endtime.julday}"
            if self.fmin or self.fmax:
                _fmin = self.fmin or 0
                _fmax = self.fmax or self.st[0].stats.sampling_rate / 2
                title += f" [{_fmin}, {_fmax}]Hz"
            # Append some information on the TauP arrivals
            if self.tp_phases:
                title += (f"\n(TauP={self.tp_model}; $\\Delta$="
                          f"{self.tp_dist_deg}deg; Z={self.tp_depth}km)")
            if self.title_append:
                title += f"\n{self.title_append}"
        else:
            title = self.title
        plt.title(title)
        # plt.suptitle(title)

        # Brute force turn off everything
        if self.frameoff:
            self.ax.set_axis_off()
            self.ax.get_legend().remove()
            self.ax.set_title("")

        plt.tight_layout()

    def finalize(self, transparent=True):
        """
        Final plot adjustments
        """
        if self.save:
            if self.save == "auto":
                _fid_out = f"{self.fids[0]}.png"
            else:
                _fid_out = self.save
            print(f"\tsaving to {_fid_out}")
            plt.savefig(_fid_out, transparent=transparent)

        if self.show:
            plt.show()
        plt.close("all")

        if self.output:
            if not os.path.exists(self.output):
                os.makedirs(self.output)
            self.st.write(os.path.join(self.output, self.fids), format="MSEED")

    def main(self):
        """
        Main plotting function
        """
        self.setup_plot()
        self.trim_waveform()
        self.process_waveforms()
        self.plot_waveforms()
        self.plot_taup_arrivals()
        if self.add_trace:
            self.plot_additional_traces()
        if self.spectrogram:
            self.plot_spectrogram()
        if self.tmarks:
            self.plot_tmarks()
        self.plot_aesthetics()
        self.finalize()


if __name__ == "__main__":
    # Initialize class using command line arguments 
    args = parse_args()
    args_dict = vars(args).copy()
    fids = args_dict.pop("fids")
    args_dict["legend"] = not args_dict.pop("no_legend", False)
    args_dict["show"] = not args_dict.pop("noshow", False)
    args_dict["sp_dbscale"] = not args_dict.pop("sp_dbscale_off", False)
    pp = PrettyPlot(fids, **args_dict)
    pp.main()
