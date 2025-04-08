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
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.dates import date2num
from obspy import read, UTCDateTime
from obspy.taup import TauPyModel
from obspy.geodetics import kilometers2degrees
from obspy.imaging.util import _set_xaxis_obspy_dates

from pysep import read_sem


SECONDS_PER_DAY = 3600.0 * 24.0


def parse_args():
    """All modifications are accomplished with command line arguments"""
    parser = argparse.ArgumentParser()

    # Waveform Processing
    parser.add_argument("fid", nargs="?", help="required, file ID")
    parser.add_argument("-tp", "--taper", nargs="?", type=float, default=0,
                        help="optional taper percentange")
    parser.add_argument("-f1", "--fmin", nargs="?", type=float, default=None,
                        help="optional filtering freqmin in Hz")
    parser.add_argument("-f2", "--fmax", nargs="?", type=float, default=None,
                        help="optional filtering freqmax in Hz")
    parser.add_argument("-z", "--zerophase", action="store_true", default=False,
                        help="apply zerophase filter or not")
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

    # TauP Phase Arrivals
    parser.add_argument("--tp_phases", nargs="+", type=str, default=None,
                        help="taup phase arrivals, requires all 'tp_*' params")
    parser.add_argument("--tp_model", nargs="?", type=str, default="iasp91",
                        help="taup model, defaults to 'iasp91'")
    parser.add_argument("--tp_dist", nargs="?", type=float, default=None,
                        help="TauP source receiver distance in km")
    parser.add_argument("--tp_depth", nargs="?", type=float, default=None,
                        help="TauP source depth km")

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
    parser.add_argument("-c", "--color", nargs="?", type=str, default="k",
                        help="color of the time series line")
    parser.add_argument("-lw", "--linewidth", nargs="?", type=float, default=0.5,
                        help="linewidth of the time series line")
    parser.add_argument("--ylabel", nargs="?", type=str, default=None,
                        help="label for units, defaults to displacement")
    parser.add_argument("--title", nargs="?", type=str, default=None,
                        help="title of the figure, defaults to ID and fmin/max")
    parser.add_argument("-ta", "--title_append", nargs="?", type=str, 
                        default="", help="append to default title")
    parser.add_argument("-tm", "--tmarks", nargs="+", 
                        help="plot vertical lines at given relative times, "
                             "should match the units of `time`. If `time`=='a' "
                             "then each tmark should be a datetime "
                             "YYYY-MM-DDTHH:MM:SS")

    # Misc
    parser.add_argument("-s", "--save", type=str, default=None,
                        help="filename to save figure")
    parser.add_argument("--noshow", action="store_true", default=False,
                        help="dont show the figure, default behavior will show")

    return parser.parse_args()

def set_plot_aesthetic(
        ax, ytick_fontsize=9., xtick_fontsize=9., tick_linewidth=1.5,
        tick_length=5., tick_direction="in", xlabel_fontsize=10.,
        ylabel_fontsize=10., axis_linewidth=1.5, spine_zorder=8, spine_top=True,
        spine_bot=True, spine_left=True, spine_right=True, title_fontsize=10.,
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
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

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


if __name__ == "__main__":
    args = parse_args()

    if not args.fid:
        sys.exit("positional argument `fid` required")

    try:
        st = read(args.fid)
    except TypeError:
        st = read_sem(args.fid)

    # Get phase arrivals from TauP if requested
    arrivals = None
    if args.tp_phases:
        assert(args.tp_dist is not None)
        assert(args.tp_depth is not None)
        dist_deg = kilometers2degrees(args.tp_dist)
        model = TauPyModel(model="iasp91")
        tp_arrivals = model.get_travel_times(source_depth_in_km=args.tp_depth,
                                             distance_in_degree=dist_deg,
                                             phase_list=args.tp_phases)

        # If some arrivals have multiple entires, only take first and last to 
        # get a range which we will plot with a window
        arrivals = {arrival.name: [] for arrival in tp_arrivals}
        for i, arrival in enumerate(tp_arrivals):
            arrivals[arrival.name].append(arrival.time)

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
                  zerophase=args.zerophase)
    elif args.fmin and not args.fmax:
        print(f"highpass {args.fmin}hz")
        st.filter("highpass", freq=args.fmin, zerophase=args.zerophase)
    elif args.fmax and args.fmin is None:
        print(f"lowpass {args.fmax}hz")
        st.filter("lowpass", freq=args.fmax, zerophase=args.zerophase)
    if args.fmin or args.fmax:
        print(f"zerophase={args.zerophase}")

    # Main plotting start
    f, ax = plt.subplots(figsize=(8, 4), dpi=200)

    if args.time.startswith("a"):
        # convert seconds of relative sample times to days and add
        # start time of trace.
        # Allow time shift a+i or a-i where i is time in hours
        sign = ""; shift = ""
        if len(args.time) > 1:
            sign = args.time[1]
            shift = args.time[2:]
            int_sign = {"+": 1, "-": -1}[sign]
            shift = int(shift)
            st[0].stats.starttime += int_sign * 60 * 60 * shift
        xvals = ((st[0].times() / SECONDS_PER_DAY) +
                        date2num(st[0].stats.starttime.datetime))
        _set_xaxis_obspy_dates(ax)
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


    plt.plot(xvals, st[0].data, c=args.color, lw=args.linewidth, zorder=6)

    # Plot phases from TauP
    if arrivals:
        for i, (name, times) in enumerate(arrivals.items()):
            if times[0] == times[-1]:
                alpha = 1
            else: 
                alpha = 0.3
            plt.axvspan(times[0], times[-1], label=name, color=f"C{i}", 
                        alpha=alpha, zorder=7)
        plt.legend(fontsize=8, loc="upper left", frameon=False)

    # Set plot aesthetics
    if args.time.startswith("a"):
        plt.xlabel(f"Time [UTC{sign}{shift}]")
    else:
        plt.xlabel(f"Time [{args.time}]")
    plt.ylabel(args.ylabel or "Displacement [m]")

    # Subset x axis
    xstart, xend = xvals.min(), xvals.max()
    if args.xlim:
        if args.time == "a":
            xstart = date2num(UTCDateTime(args.xlim[0]).datetime)
            xend = date2num(UTCDateTime(args.xlim[1]).datetime)
        else:
            xstart, xend = args.xlim
    plt.xlim(xstart, xend)

    # Subset y axis
    ymax = np.amax([st[0].data.min(), st[0].data.max()])
    ymin = -1 * ymax
    if args.ylim:
        # Allow for one entry to set min/max if they're the same
        if len(args.ylim) == 1:
            ymin, ymax = [-1 * args.ylim[0], args.ylim[0]]
        else:
            ymin, ymax = args.ylim
    else:
        # Need to rescale y-axis based on the subset x-axis, hacky
        idx_start = (np.abs(xvals - xstart)).argmin()
        idx_end = (np.abs(xvals - xend)).argmin()
        yvals = st[0].data[idx_start: idx_end]
        ymax = np.amax([yvals.min(), yvals.max()])
        ymin = -1 * ymax

    plt.ylim(ymin, ymax)

    # Add vertical lies at certain times
    if args.tmarks:
        for tmark in args.tmarks:
            if args.time == "a":
                tmark = date2num(UTCDateTime(tmark).datetime)
                # date2num(st[0].stats.starttime.datetime) + tmark
            plt.axvline(tmark, c="r", lw=0.5)

    # Finish off by setting plot aesthetics
    if not args.title:
        title = f"{st[0].get_id()} [{args.fmin}, {args.fmax}]Hz"

        # Append some information on the TauP arrivals
        if arrivals:
            title += (f"\n(TauP={args.tp_model}; $\\Delta$={args.tp_dist}km; "
                      f"Z={args.tp_depth}km)")

        title += f"\n{args.title_append}"

    else:
        title = args.title
    plt.title(title)

    set_plot_aesthetic(ax)
    f.tight_layout()

    # Finalize Plot
    if args.save == "auto":
        plt.savefig(f"{args.fid}.png")
    elif args.save is not None:
        plt.savefig(args.save)

    if not args.noshow:
        plt.show()
    plt.close("all")

