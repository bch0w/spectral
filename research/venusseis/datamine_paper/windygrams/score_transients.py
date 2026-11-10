#!/usr/bin/env python
"""
Rank 4 h wind-state windows by the impulsive transients they contain.

Runs a recursive STA/LTA on the vertical component of all co-located
sensors of a site (Covered, Uncovered, Windscreen, Permanent) and clusters
triggers across sensors. A cluster counts as a ground transient if the
permanent station and at least one of the Covered and Windscreen nodes
triggered. For the detection rate of sensor X, the requirement is instead
the permanent station plus another node besides X, so each sensor is
measured against independent stations.

Writes to `--path_out`:

- ranked_windows.csv: one row per window, sorted by site, wind state and
  number of transients. `events` = ground transients; `n_X`, `hit_X` =
  transients and how many sensor X detected; `only_X` = triggers on X
  alone (mostly wind); `cft_P`, `cft_U` = median peak STA/LTA at transients;
  `vel_P`, `vel_U` = largest 2-20 Hz velocity at transients (m/s); `rho` =
  Spearman correlation of the Uncovered 1-min amplitude and wind
- event_times.csv: onset of each ground transient, as `offset_s` seconds
  after the window start. Read by plot_wind_states.py to mark transients

Run through the symlink in the windygrams working directory. The default
output directory is linked to the copies of both files in the repo.

Example::

    python score_transients.py
    python score_transients.py --sites MITD --path_out CHECKS/test
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd
from obspy import read, read_inventory
from obspy.signal.trigger import recursive_sta_lta, trigger_onset
from scipy.stats import spearmanr

import find_windows as fw
from score_states import HOURS, PERM, STATES, candidates, minute_wind

PATH_OUT = "CHECKS/transients"
NCAND = 15
PAD = 600.
BAND = (2., 20.)
STA, LTA, ON, OFF = 1., 30., 5., 1.5
TOL = 3.  # seconds, co-located sensors have negligible moveout


def station_cft(code, t0, t1, inv_xx):
    """
    Return the characteristic function and triggers of one station

    :type code: str
    :param code: trace ID
    :type t0: obspy.UTCDateTime
    :param t0: window start
    :type t1: obspy.UTCDateTime
    :param t1: window end
    :type inv_xx: obspy.Inventory
    :param inv_xx: response information of the nodes
    :rtype: tuple
    :return: times relative to t0, STA/LTA, velocity and trigger on times
    """
    fid = f"DATA/{code}.{t0.year}.{t0.julday:03d}"
    st = read(fid, starttime=t0 - PAD, endtime=t1 + PAD)
    for tr in st:
        tr.data = tr.data.astype(float)
    tr = st.merge(fill_value="interpolate").select(
        channel=code.split(".")[-1])[0]
    fs = tr.stats.sampling_rate
    tr.detrend("linear").taper(max_percentage=None, max_length=30.)
    tr.filter("bandpass", freqmin=BAND[0], freqmax=BAND[1], corners=4)
    cft = recursive_sta_lta(tr.data, int(STA * fs), int(LTA * fs))
    times = (tr.stats.starttime - t0) + np.arange(tr.stats.npts) / fs
    onsets = np.array([times[i] for i, _ in trigger_onset(cft, ON, OFF)])
    if len(onsets):
        onsets = onsets[(onsets >= 0) & (onsets <= t1 - t0)]
    if code.startswith("XX."):
        inv = inv_xx
    else:
        inv = read_inventory(f"XML/{'.'.join(code.split('.')[:2])}.xml")
    vel = tr.data / inv.get_response(code, t0).instrument_sensitivity.value
    return times, cft, vel, onsets


def clusters(onsets):
    """
    Group trigger times of all stations into clusters

    :type onsets: dict
    :param onsets: trigger on times keyed by station letter
    :rtype: list of tuple
    :return: (first time, set of stations) of each cluster
    """
    triggers = sorted((t, sta) for sta, on in onsets.items() for t in on)
    out = []
    for t, sta in triggers:
        if out and t - out[-1][2] <= TOL:
            out[-1][1].add(sta)
            out[-1][2] = t
        else:
            out.append([t, {sta}, t])
    return [(first, stas) for first, stas, _ in out]


def window_peak(times, values, center, half=5.):
    """
    Return the largest absolute value of a series around a time

    :type times: numpy.ndarray
    :param times: time axis
    :type values: numpy.ndarray
    :param values: series
    :type center: float
    :param center: center time
    :type half: float
    :param half: half width in seconds
    :rtype: float
    :return: peak value
    """
    i0, i1 = np.searchsorted(times, [center - half, center + half])
    return np.abs(values[i0:i1]).max() if i1 > i0 else np.nan


def wind_correlation(res, site, t0, t1):
    """
    Return the correlation of the Uncovered node's amplitude and the wind

    :type res: tuple
    :param res: output of `station_cft` for the Uncovered node
    :type site: str
    :param site: site name
    :type t0: obspy.UTCDateTime
    :param t0: window start
    :type t1: obspy.UTCDateTime
    :param t1: window end
    :rtype: float
    :return: Spearman correlation of 1-min mean amplitude and wind speed
    """
    times, _, vel, _ = res
    keep = (times >= 0) & (times <= t1 - t0)
    nmin = int((t1 - t0) // 60)
    bins = (times[keep] // 60).astype(int).clip(0, nmin - 1)
    total = np.bincount(bins, np.abs(vel[keep]), nmin)
    env = total / np.maximum(np.bincount(bins, minlength=nmin), 1)
    return spearmanr(env, minute_wind(site, t0, t1, nmin))[0]


def score_window(site, cha, t0, t1, inv_xx):
    """
    Detect the transients of a window on all co-located sensors

    :type site: str
    :param site: site name
    :type cha: str
    :param cha: node channel code
    :type t0: obspy.UTCDateTime
    :param t0: window start
    :type t1: obspy.UTCDateTime
    :param t1: window end
    :type inv_xx: obspy.Inventory
    :param inv_xx: response information of the nodes
    :rtype: dict
    :return: detection statistics and `cl_times`, the transient onsets
    """
    codes = {sta: f"XX.{site}{sta}..{cha}" for sta in "CUW"}
    codes["P"] = PERM[site]
    res = {sta: station_cft(code, t0, t1, inv_xx)
           for sta, code in codes.items()}
    groups = clusters({sta: val[3] for sta, val in res.items()})
    row = {}
    for sta in "CUW":
        others = set("CUW") - {sta}
        ground = [g for g in groups if "P" in g[1] and g[1] & others]
        row[f"n_{sta}"] = len(ground)
        row[f"hit_{sta}"] = sum(sta in g[1] for g in ground)
        row[f"only_{sta}"] = sum(g[1] == {sta} for g in groups)
    # Independent of the Uncovered node, these are used for the ranking
    ground = [g for g in groups if "P" in g[1] and g[1] & set("CW")]
    row["events"] = len(ground)
    row["cl_times"] = [round(g[0]) for g in ground]
    if ground:
        for name, col, func in (("cft", 1, np.median), ("vel", 2, np.max)):
            for sta in "PU":
                row[f"{name}_{sta}"] = func([
                    window_peak(res[sta][0], res[sta][col], g[0])
                    for g in ground])
    row["rho"] = wind_correlation(res["U"], site, t0, t1)
    return row


def write_outputs(rows, path_out):
    """
    Write the ranked windows and transient onset times

    :type rows: list of dict
    :param rows: statistics of every window
    :type path_out: str
    :param path_out: output directory
    """
    os.makedirs(path_out, exist_ok=True)
    df = pd.DataFrame(rows)
    events = pd.DataFrame([
        {"site": row.site, "jday": row.jday, "start": row.t0,
         "state": row.state, "offset_s": offset}
        for row in df.itertuples() for offset in row.cl_times])
    events.drop_duplicates(["site", "jday", "start", "offset_s"]).to_csv(
        os.path.join(path_out, "event_times.csv"), index=False)
    ranked = df.drop(columns="cl_times").sort_values(
        ["site", "state", "events", "hit_U"],
        ascending=[True, True, False, False])
    ranked.to_csv(os.path.join(path_out, "ranked_windows.csv"), index=False,
                  float_format="%.4g")


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
    """Score the candidate windows of every site and write the outputs"""
    args = parse_args()
    warnings.filterwarnings("ignore")
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    inv_xx = read_inventory("XML/ISHTAR.xml")
    rows = []
    for site, info in fw.SITES.items():
        if args.sites and site not in args.sites:
            continue
        df = fw.rank_windows(fw.wind_10min(site), HOURS)
        for state, key, ascending in STATES:
            for _, row, t0, t1 in candidates(site, df, key, ascending, NCAND):
                rows.append({
                    "site": site, "state": state, "jday": t0.julday,
                    "t0": str(t0)[11:16], "t1": str(t1)[11:16],
                    "wmean": row["mean"], "wmin": row["min"],
                    "wmax": row["max"], "d": row["delta"],
                    **score_window(site, info["cha"], t0, t1, inv_xx)})
                print(site, state, rows[-1]["jday"], rows[-1]["t0"],
                      "events", rows[-1]["events"], flush=True)
    write_outputs(rows, args.path_out)


if __name__ == "__main__":
    main()
