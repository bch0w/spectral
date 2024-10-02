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
from prettyplot import set_plot_aesthetic

import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy as np

from obspy import read
from pysep import read_sem
from scipy.fftpack import fft, fftfreq, next_fast_len

def parse_args():
    """All modifications are accomplished with command line arguments"""
    parser = argparse.ArgumentParser()

    # Waveform Processing
    parser.add_argument("fid", nargs="?", help="required, file ID")

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

    # Separate figures
    f, ax = plt.subplots(figsize=(4, 4), dpi=200)

    # Process waveform for spectra
    tr = st[0]  # assuming only ONE trace in the stream
    tr.detrend("demean")
    sr = tr.stats.sampling_rate
    nfft = next_fast_len(tr.stats.npts)  # pad data with zeros for fast FT
    pos_freq = (nfft + 1) // 2  # index for positive frequencies
    spec = fft(tr.data, nfft)[:pos_freq]
    freq = fftfreq(nfft, 1 / sr)[:pos_freq]  # get freqs of DFT bin


    plt.plot(freq, np.abs(spec), c=args.color, lw=args.linewidth)  


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

