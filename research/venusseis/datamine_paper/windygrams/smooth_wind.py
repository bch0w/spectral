#!/usr/bin/env python
"""
Write running-mean wind speed files for plotting with ``pplot --add_trace``.

The 1 Hz anemometer speed (LWS) fills the waveform panel as a solid band
over a few hours, so it is averaged with a centred running mean and saved
as MiniSEED next to the raw data, e.g.
``WIND/XX.K24KA..LWS.2025.267`` from ``DATA/XX.K24KA..LWS.2025.267``.

Example::

    python smooth_wind.py DATA/XX.K24KA..LWS.2025.267 --window 60
"""
import argparse
import os

import numpy as np
from obspy import read


def smooth(fid, window=60.0, path_out="WIND"):
    """
    Merge one wind speed file, smooth it with a running mean and save it

    :type fid: str
    :param fid: path to an ObsPy-readable wind speed file
    :type window: float
    :param window: running mean length in seconds
    :type path_out: str
    :param path_out: directory to write the smoothed file to
    :rtype: str
    :return: path of the smoothed file
    """
    st = read(fid).merge(method=1, fill_value="interpolate")
    tr = st[0]
    n = max(int(round(window * tr.stats.sampling_rate)), 1)
    data = tr.data.astype(np.float64)
    tr.data = np.convolve(data, np.ones(n) / n, mode="same")
    # Edges average over fewer samples; renormalize so they are not biased
    tr.data /= np.convolve(np.ones_like(data), np.ones(n) / n, mode="same")
    os.makedirs(path_out, exist_ok=True)
    fid_out = os.path.join(path_out, os.path.basename(fid))
    tr.write(fid_out, format="MSEED")
    return fid_out


def main():
    """Parse command line arguments and smooth each file given"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fids", nargs="+", help="wind speed file(s)")
    parser.add_argument("--window", type=float, default=60.0,
                        help="running mean length in seconds, default 60")
    parser.add_argument("--path_out", default="WIND",
                        help="output directory, default WIND")
    args = parser.parse_args()
    for fid in args.fids:
        print(smooth(fid, args.window, args.path_out))


if __name__ == "__main__":
    main()
