"""
Plot spectrogram alongside waveform
"""
from prettyplot import set_plot_aesthetic

import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy as np

from obspy import read, UTCDateTime
from pysep import read_sem
from scipy.fftpack import fft, fftfreq, next_fast_len

def parse_args():
    """All modifications are accomplished with command line arguments"""
    parser = argparse.ArgumentParser()

    # Waveform Processing
    parser.add_argument("fid", nargs="?", help="required, file ID")

    parser.add_argument("-t0", "--tstart", nargs="?", type=float, default=None,
                        help="start time in UTCDateTime format to cut waveform")
    parser.add_argument("-t1", "--tend", nargs="?", type=float, default=None,
                        help="end time in UTCDateTime to cut waveform")

    # Plot Aesthetics
    parser.add_argument("-x", "--xlim", nargs="+", type=float, default=None,
                        help="time axis limits in s")
    parser.add_argument("-y", "--ylim", nargs="+", type=float, default=None,
                        help="time axis limits in s")
    parser.add_argument("-c", "--color", nargs="?", type=str, default="k",
                        help="color of the time series line")
    parser.add_argument("-lw", "--linewidth", nargs="?", type=float, default=0.5,
                        help="linewidth of the time series line")
    parser.add_argument("--title", nargs="?", type=str, default=None,
                        help="title of the figure, defaults to ID and fmin/max")
    parser.add_argument("-ta", "--title_append", nargs="?", type=str, 
                        default="", help="append to default title")

    # Misc
    parser.add_argument("-s", "--save", type=str, default=None,
                        help="filename to save figure")
    parser.add_argument("-S", "--spectra", action="store_true", default=False,
                        help="plot spectra of the raw trace")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if not args.fid:
        sys.exit("positional argument `fid` required")

    try:
        st = read(args.fid)
    except TypeError:
        st = read_sem(args.fid)

    # Trim waveform
    if args.tstart:
        tstart = UTCDateTime(args.tstart)
    else:
        tstart = args.tstart
    if args.tend:
        tend = UTCDateTime(args.tend)
    else:
        tend = args.tend
    st.trim(tstart, tend)

    st.spectrogram()
    



    # Set plot aesthetics
    plt.xlabel("Freq. [Hz]")
    plt.ylabel("Amplitude")

    if args.xlim:
        plt.xlim(args.xlim)

    if args.ylim:
        plt.ylim(args.ylim)

    if not args.title:
        title = f"{st[0].get_id()} Amplitude Spectra"
        title += f"\n{args.title_append}"
    else:
        title = args.title

    plt.title(title)

    set_plot_aesthetic(ax)
    f.tight_layout()

    # Finalize Plot
    if args.save:
        plt.savefig(args.save)

    plt.show()
    plt.close("all")

