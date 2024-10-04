"""
Hijack ObsPy Dayplot and make it look a bit nicer

.. rubric

    python dayplot.py <fid> --fmin 1 --fmax 100 
"""
import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from obspy import read


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
    parser.add_argument("--integrate", nargs="?", type=int, default=0,
                        help="Integrate the time series, will demean and taper,"
                             " value for integrate will be the number of times")
    parser.add_argument("--differentiate", nargs="?", type=int, default=0,
                        help="Diff. the time series, will demean and taper,"
                             " value for integrate will be the number of times")
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
    # Thicken up the bounding axis lines
    for axis, flag in zip(["top", "bottom", "left", "right"],
                          [spine_top, spine_bot, spine_left, spine_right]):
        # Deal with the case where command line users are inputting strings
        if isinstance(flag, str):
            flag = bool(flag.capitalize() == "True")
        ax.spines[axis].set_visible(flag)
        ax.spines[axis].set_linewidth(axis_linewidth)

    # Remove axis labels these will be replaced by times
    ax.xaxis.label.set_size(0)
    ax.yaxis.label.set_size(0)
    ax.tick_params(axis="x", labelsize=0)
    ax.tick_params(axis="y", labelsize=0)


if __name__ == "__main__":
    args = parse_args()

    try:
        st = read(args.fid)
    except TypeError:
        st = read_sem(args.fid)


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

    f, ax = plt.subplots(figsize=(6, 8), dpi=100)
    set_plot_aesthetic(ax)

    st.plot(type="dayplot", interval=60, right_vertical_labels=False,
            linewidth=0.25, one_tick_per_line=True, method="full",
            color=[f"C{_}" for _ in range(0, 6)], show_y_UTC_label=False, 
            show=False, fig=f, title=os.path.basename(args.fid))
    f.tight_layout()

    if args.save:
        plt.savefig(args.save)

    if not args.noshow:
        plt.show()



