#!/usr/bin/env python
"""
Score candidate 4 h wind-state windows at every site.

Takes the highest ranked windows of each wind state from find_windows.py
(high, low, rising, falling; at least 12 h apart, full data coverage at the
nodes and the permanent station) and scores how well the 5 Hz highpassed
amplitude of the Uncovered node follows the wind. Used to choose the panels
of the first wind-state figures (`STATES` in plot_wind_states.py).

Columns of the output CSV: wind statistics of 10-min means (`wmean`, `wmin`,
`wmax`, `d` = last minus first hour); median and peak 1-min amplitude
(`envmed`, `envmax`, m/s); `ratio` = median amplitude of the last hour over
the first hour; `rho` = Spearman correlation of 1-min amplitude and wind;
`tran` = fraction of minutes louder than 3x the running median amplitude;
`lst` = local start time.

Run through the symlink in the windygrams working directory.

Example::

    python score_states.py
    python score_states.py --sites MITD K24K
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd
from obspy import UTCDateTime, read, read_inventory
from scipy.fft import next_fast_len
from scipy.ndimage import median_filter, uniform_filter1d
from scipy.signal import hilbert
from scipy.stats import spearmanr

import find_windows as fw

PATH_OUT = "CHECKS/wind_states"
PERM = {"DESD": "HV.DESD..EHZ", "H23K": "AK.H23K..BHZ", "K24K": "AK.K24K..BHZ",
        "KAED": "HV.KAED..EHZ", "MDM": "AK.MDM..BHZ", "MITD": "HV.MITD..HHZ",
        "PAX": "AK.PAX..BHZ", "RND": "AK.RND..BHZ"}
HOURS = 4
NCAND = 6
MIN_SEP = 12 * 3600
# (wind state, find_windows column ranked, ascending order)
STATES = (("high", "min", False), ("low", "max", True),
          ("rising", "delta", False), ("falling", "delta", True))


def has_permanent(site, t0, t1):
    """
    Check that the permanent station of a site covers a window

    :type site: str
    :param site: site name
    :type t0: obspy.UTCDateTime
    :param t0: window start
    :type t1: obspy.UTCDateTime
    :param t1: window end
    :rtype: bool
    :return: True if gaps inside the window total less than 30 s
    """
    fid = f"DATA/{PERM[site]}.{t0.year}.{t0.julday:03d}"
    spans = fw.gaps_in_day(fid)
    cov = sum(max(0., min(e, t1) - max(s, t0)) for s, e in spans or [])
    return cov >= (t1 - t0) - 30


def candidates(site, df, key, ascending, ncand):
    """
    Yield the best windows of one ranking that pass the coverage checks

    :type site: str
    :param site: site name
    :type df: pandas.DataFrame
    :param df: output of find_windows.rank_windows
    :type key: str
    :param key: column to rank by
    :type ascending: bool
    :param ascending: rank in ascending order
    :type ncand: int
    :param ncand: largest number of windows to yield
    :rtype: tuple
    :return: window start time, its row in `df`, UTCDateTime start and end
    """
    used = []
    for t, row in df.sort_values(key, ascending=ascending).iterrows():
        if any(abs((t - u).total_seconds()) < MIN_SEP for u in used):
            continue
        t0 = UTCDateTime(t.to_pydatetime())
        t1 = t0 + HOURS * 3600
        if not fw.coverage_ok(site, t0, t1, max_gap=30.):
            continue
        if not has_permanent(site, t0, t1):
            continue
        used.append(t)
        yield t, row, t0, t1
        if len(used) == ncand:
            return


def minute_wind(site, t0, t1, nmin):
    """
    Return the mean wind speed of each minute of a window

    :type site: str
    :param site: site name
    :type t0: obspy.UTCDateTime
    :param t0: window start
    :type t1: obspy.UTCDateTime
    :param t1: window end
    :type nmin: int
    :param nmin: number of minutes
    :rtype: numpy.ndarray
    :return: 1-min mean wind speed in m/s
    """
    st = read(f"DATA/XX.{site}A..LWS.{t0.year}.{t0.julday:03d}",
              starttime=t0, endtime=t1)
    for tr in st:
        tr.data = tr.data.astype(float)  # some PAX days have mixed dtypes
    wind = st.merge(fill_value="interpolate")[0]
    return np.array([wind.data[int(i * 60):int((i + 1) * 60)].mean()
                     for i in range(nmin)])


def score_window(site, cha, t0, t1, inv):
    """
    Score how well the Uncovered node's amplitude follows the wind

    :type site: str
    :param site: site name
    :type cha: str
    :param cha: node channel code
    :type t0: obspy.UTCDateTime
    :param t0: window start
    :type t1: obspy.UTCDateTime
    :param t1: window end
    :type inv: obspy.Inventory
    :param inv: node response information
    :rtype: dict
    :return: 'envmed', 'envmax', 'ratio', 'rho' and 'tran'
    """
    code = f"XX.{site}U..{cha}"
    tr = read(f"DATA/{code}.{t0.year}.{t0.julday:03d}",
              starttime=t0 - 300, endtime=t1 + 300).merge(
                  fill_value="interpolate").select(channel=cha)[0]
    tr.detrend("linear")
    tr.filter("highpass", freq=5, corners=4)
    data = tr.data / inv.get_response(code, t0).instrument_sensitivity.value
    env = np.abs(hilbert(data, N=next_fast_len(len(data))))[:len(data)]
    n60 = int(60 * tr.stats.sampling_rate)
    env = uniform_filter1d(env, n60)[::n60]
    times = -300 + np.arange(len(env)) * 60 + 30.
    env = env[(times >= 0) & (times <= t1 - t0)]
    wind = minute_wind(site, t0, t1, len(env))
    npts = min(len(wind), len(env))
    base = median_filter(env, size=31, mode="nearest")
    return {"envmed": np.median(env), "envmax": env.max(),
            "ratio": np.median(env[-60:]) / np.median(env[:60]),
            "rho": spearmanr(env[:npts], wind[:npts])[0],
            "tran": np.mean(env > 3 * base)}


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites", nargs="+", default=None,
                        help="sites to score, default all")
    parser.add_argument("--path_out", default=PATH_OUT,
                        help="output directory")
    return parser.parse_args()


def main():
    """Score the candidate windows of every site and write a CSV"""
    args = parse_args()
    warnings.filterwarnings("ignore")
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    inv = read_inventory("XML/ISHTAR.xml")
    rows = []
    for site, info in fw.SITES.items():
        if args.sites and site not in args.sites:
            continue
        df = fw.rank_windows(fw.wind_10min(site), HOURS)
        for state, key, ascending in STATES:
            for t, row, t0, t1 in candidates(site, df, key, ascending, NCAND):
                rows.append({
                    "site": site, "state": state, "jday": t0.julday,
                    "t0": str(t0)[11:16], "t1": str(t1)[11:16],
                    "wmean": row["mean"], "wmin": row["min"],
                    "wmax": row["max"], "d": row["delta"],
                    **score_window(site, info["cha"], t0, t1, inv),
                    "lst": t.tz_convert(info["tz"]).strftime("%H:%M")})
    os.makedirs(args.path_out, exist_ok=True)
    fid = os.path.join(args.path_out, "scored_windows.csv")
    pd.DataFrame(rows).to_csv(fid, index=False, float_format="%.6g")
    print(fid)


if __name__ == "__main__":
    main()
