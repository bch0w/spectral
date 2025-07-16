"""
Plot spectrogram alongside waveform
"""
from prettyplot import set_plot_aesthetic

import argparse
import sys
import os
import math
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import mlab
from matplotlib.colors import Normalize
from obspy import read, UTCDateTime
from obspy.imaging.cm import obspy_sequential
from pysep import read_sem
from scipy.fftpack import fft, fftfreq, next_fast_len


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


def spectrogram(data, samp_rate, per_lap=0.9, wlen=None, log=False,
                outfile=None, fmt=None, axes=None, dbscale=False,
                mult=8.0, cmap=obspy_sequential, zorder=None, title=None,
                show=True, clip=[0.0, 1.0]):
    """
    Modified from ObsPys Spectrogram function
    https://docs.obspy.org/_modules/obspy/imaging/spectrogram.html#spectrogram

    Computes and plots spectrogram of the input data.

    :param data: Input data
    :type samp_rate: float
    :param samp_rate: Samplerate in Hz
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
    # enforce float for samp_rate
    samp_rate = float(samp_rate)

    # set wlen from samp_rate if not specified otherwise
    if not wlen:
        wlen = 128 / samp_rate

    npts = len(data)

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

    data = data - data.mean()
    end = npts / samp_rate

    # Here we call not plt.specgram as this already produces a plot
    # matplotlib.mlab.specgram should be faster as it computes only the
    # arrays
    # XXX mlab.specgram uses fft, would be better and faster use rfft
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

    if not axes:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        ax = axes

    # calculate half bin width
    halfbin_time = (time[1] - time[0]) / 2.0
    halfbin_freq = (freq[1] - freq[0]) / 2.0

    kwargs = {'cmap': cmap, 'zorder': zorder}
    if log:
        # pcolor expects one bin more at the right end
        freq = np.concatenate((freq, [freq[-1] + 2 * halfbin_freq]))
        time = np.concatenate((time, [time[-1] + 2 * halfbin_time]))
        # center bin
        time -= halfbin_time
        freq -= halfbin_freq
        # Log scaling for frequency values (y-axis)
        ax.set_yscale('log')
        # Plot times
        ax.pcolormesh(time, freq, specgram, norm=norm, **kwargs)
    else:
        # this method is much much faster!
        specgram = np.flipud(specgram)
        # center bin
        extent = (time[0] - halfbin_time, time[-1] + halfbin_time,
                  freq[0] - halfbin_freq, freq[-1] + halfbin_freq)
        ax.imshow(specgram, interpolation="nearest", extent=extent, **kwargs)

    # set correct way of axis, whitespace before and after with window
    # length
    ax.axis('tight')
    ax.set_xlim(0, end)
    ax.grid(False)

    if axes:
        return ax

    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Frequency [Hz]')
    if title:
        ax.set_title(title)

    return fig

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

    spectrogram(st[0].data, st[0].stats.sampling_rate, log=True) 
    plt.show()



    a=1/0
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

