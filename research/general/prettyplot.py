"""
Pretty plotting of ObsPy Streams
"""
import argparse
import sys
import os
import matplotlib.pyplot as plt

from obspy import read
from pysep import read_sem

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("fid", nargs=1)
    parser.add_argument("--xlim", nargs="+", type=int, default=None)
    parser.add_argument("--taper", nargs="?", type=float, default=0)
    parser.add_argument("--fmin", nargs="?", type=int, default=None)
    parser.add_argument("--fmax", nargs="?", type=int, default=None)
    parser.add_argument("--t0", nargs="?", type=float, default=0)

    parser.add_argument("--c", nargs="?", type=str, default="k")
    parser.add_argument("--lw", nargs="?", type=float, default=1)

    return parser

def set_plot_aesthetic(
        ax, ytick_fontsize=8., xtick_fontsize=12., tick_linewidth=1.5,
        tick_length=5., tick_direction="in", xlabel_fontsize=10.,
        ylabel_fontsize=10., axis_linewidth=2., spine_zorder=8, spine_top=True,
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
        plt.grid(visible=True, which="major", axis="x", alpha=0.5, linewidth=1)
    if xgrid_minor:
        plt.grid(visible=True, which="minor", axis="x", alpha=0.2, linewidth=.5)
    if ygrid_major:
        plt.grid(visible=True, which="major", axis="y", alpha=0.5, linewidth=1)
    if ygrid_minor:
        plt.grid(visible=True, which="minor", axis="y", alpha=0.2, linewidth=.5)


if __name__ == "__main__":
    args = parse_args()

    try:
        st = read(args.fid)
    except TypeError:
        st = read_sem(args.fid)

    f, ax = plt.subplot()

    if args.taper:
        st.taper(args.taper)
    if fmin:
        st.filter("bandpass", freqmin=fmin, freqmax=fmax)

    plt.plot(st[0].times - args.t0, st[0].data, c=args.c, lw=args.lw)
    if args.xlim:
        plt.set_xlim(args.xlim)

    set_plot_aesthetic(ax)
    plt.show()


