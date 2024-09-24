"""
Pretty plotting of ObsPy Streams, can either be ObsPy 'read'able format or
two-column ASCII output for SPECFEM synthetics. Most things controlled by parser
"""
import argparse
import sys
import os
import matplotlib.pyplot as plt

from obspy import read
from pysep import read_sem

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("fid", nargs="?", help="required, file ID")
    parser.add_argument("--xlim", nargs="+", type=int, default=None,
                        help="time axis limits in s")
    parser.add_argument("--taper", nargs="?", type=float, default=0,
                        help="optional taper percentange")
    parser.add_argument("--fmin", nargs="?", type=int, default=None,
                        help="optional filtering freqmin in Hz")
    parser.add_argument("--fmax", nargs="?", type=int, default=None,
                        help="optional filtering freqmax in Hz")
    parser.add_argument("--t0", nargs="?", type=float, default=0,
                        help="SPECFEM USER_T0 if synthetics")

    parser.add_argument("--c", nargs="?", type=str, default="k",
                        help="color of the time series line")
    parser.add_argument("--lw", nargs="?", type=float, default=1,
                        help="linewidth of the time series line")
    parser.add_argument("--ylabel", nargs="?", type=str, default=None,
                        help="label for units, defaults to displacement")
    parser.add_argument("--title", nargs="?", type=str, default=None,
                        help="title of the figure, defaults to ID and fmin/max")

    return parser.parse_args()

def set_plot_aesthetic(
        ax, ytick_fontsize=9., xtick_fontsize=9., tick_linewidth=1.5,
        tick_length=5., tick_direction="in", xlabel_fontsize=10.,
        ylabel_fontsize=10., axis_linewidth=2., spine_zorder=8, spine_top=True,
        spine_bot=True, spine_left=True, spine_right=True, title_fontsize=12.,
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

    try:
        st = read(args.fid)
    except TypeError:
        st = read_sem(args.fid)

    f, ax = plt.subplots(figsize=(8, 4), dpi=200)

    if args.taper:
        st.taper(args.taper)
    if args.fmin:
        st.filter("bandpass", freqmin=fmin, freqmax=fmax)

    xvals = st[0].times() - args.t0
    plt.plot(xvals, st[0].data, c=args.c, lw=args.lw)
    plt.xlabel("Time [s]")
    plt.ylabel(args.ylabel or "Displacement [m]")

    if args.xlim:
        plt.xlim(args.xlim)
    else:
        plt.xlim(xvals.min(), xvals.max())

    plt.title(args.title or f"{st[0].get_id()} [{args.fmin}, {args.fmax}]Hz")

    set_plot_aesthetic(ax)
    plt.show()


