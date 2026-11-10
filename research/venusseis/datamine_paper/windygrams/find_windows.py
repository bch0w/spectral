#!/usr/bin/env python
"""
Find high-, low- and changing-wind windows for each ISHTAR (XX) site.

Wind speed (LWS, 1 Hz) is averaged to 10-min means, then 4 h windows are
ranked as:

- high: the largest *minimum* 10-min speed, i.e. wind that stays strong
- low: the smallest *maximum* 10-min speed, i.e. calm throughout. A cup
  anemometer reads 0 below its start-up speed, so calm can read 0.0. Among
  candidates within `low_tol` m/s of the calmest, the one closest in local
  time of day to the high window is taken, so cultural noise is similar
- change: the largest difference between the first- and last-hour mean
  speed, rising or falling

A window must lie inside one UTC day (pplot reads one day file), every
seismic station (C, U, W) and the anemometer must have data for it, and
gaps inside it must total less than `max_gap` seconds. Windows are printed
as entries for the `WINDOWS` dict of windygrams.py.

Example::

    python find_windows.py                 # all sites (needs pandas)
    python find_windows.py --sites PAX RND --hours 2
"""
import argparse
import glob
import os
import re
from functools import lru_cache

import numpy as np
import pandas as pd
from obspy import read, UTCDateTime

PATH_DATA = "DATA"
SITES = {
    "DESD": {"cha": "HHZ", "tz": "Pacific/Honolulu"},
    "H23K": {"cha": "DHZ", "tz": "America/Anchorage"},
    "K24K": {"cha": "EHZ", "tz": "America/Anchorage"},
    "KAED": {"cha": "HHZ", "tz": "Pacific/Honolulu"},
    "MDM": {"cha": "DHZ", "tz": "America/Anchorage"},
    "MITD": {"cha": "HHZ", "tz": "Pacific/Honolulu"},
    "PAX": {"cha": "EHZ", "tz": "America/Anchorage"},
    "RND": {"cha": "EHZ", "tz": "America/Anchorage"},
}


def wind_10min(site, min_frac=0.5):
    """
    Read every wind speed day file of a site and average to 10-min means

    :type site: str
    :param site: site name, e.g. 'PAX'
    :type min_frac: float
    :param min_frac: minimum fraction of 1 Hz samples needed in a 10-min
        bin, otherwise it is NaN
    :rtype: pandas.Series
    :return: 10-min mean wind speed indexed by UTC bin start
    """
    parts = []
    for fid in sorted(glob.glob(f"{PATH_DATA}/XX.{site}A..LWS.*")):
        for tr in read(fid):
            t0 = pd.Timestamp(tr.stats.starttime.datetime, tz="UTC")
            idx = t0 + pd.to_timedelta(tr.times(), unit="s")
            parts.append(pd.Series(tr.data.astype(float), index=idx))
    s = pd.concat(parts).sort_index()
    s = s[~s.index.duplicated()]
    # Fault codes (e.g. 1000) come in bursts with corrupt neighbors, so a
    # bin holding any fault sample is dropped, as in plot_wind.py
    fault = (s < 0) | (s >= 100)
    grp = s.resample("10min")
    mean, count = grp.mean(), grp.count()
    mean[count < min_frac * 600] = np.nan
    mean[fault.resample("10min").sum() > 0] = np.nan
    return mean


@lru_cache(maxsize=None)
def gaps_in_day(fid):
    """
    Return the data spans of a day file, or None if it does not exist

    :type fid: str
    :param fid: path to the MiniSEED day file
    :rtype: list of (UTCDateTime, UTCDateTime) or None
    :return: start and end of each trace in the file
    """
    if not os.path.exists(fid):
        return None
    st = read(fid, headonly=True)
    return [(tr.stats.starttime, tr.stats.endtime) for tr in st]


def coverage_ok(site, start, end, max_gap=60.):
    """
    Check that all seismic stations and the anemometer cover a window

    :type site: str
    :param site: site name
    :type start: obspy.UTCDateTime
    :param start: window start
    :type end: obspy.UTCDateTime
    :param end: window end
    :type max_gap: float
    :param max_gap: largest total gap in seconds allowed in the window
    :rtype: bool
    :return: True if every channel has data for the window
    """
    jday = start.julday
    cha = SITES[site]["cha"]
    codes = [f"XX.{site}{s}..{cha}" for s in "CUW"] + [f"XX.{site}A..LWS"]
    for code in codes:
        spans = gaps_in_day(f"{PATH_DATA}/{code}.{start.year}.{jday:03d}")
        if not spans:
            return False
        covered = sum(max(0., min(e, end) - max(s, start))
                      for s, e in spans)
        if (end - start) - covered > max_gap:
            return False
    return True


def rank_windows(speed, hours=4):
    """
    Compute statistics of every window of `hours` length on the 10-min grid

    :type speed: pandas.Series
    :param speed: 10-min mean wind speed
    :type hours: int
    :param hours: window length in hours
    :rtype: pandas.DataFrame
    :return: mean, min, max and change (last minus first hour mean) per
        window, indexed by window start; windows crossing midnight UTC or
        with missing 10-min bins are dropped
    """
    n = hours * 6
    roll = speed.rolling(n, min_periods=n)
    hour = speed.rolling(6, min_periods=6).mean()
    df = pd.DataFrame({"mean": roll.mean(), "min": roll.min(),
                       "max": roll.max(),
                       "delta": hour - hour.shift(n - 6)})
    df.index = df.index - pd.Timedelta(minutes=10 * (n - 1))
    end = df.index + pd.Timedelta(hours=hours)
    same_day = (end - pd.Timedelta(seconds=1)).date == df.index.date
    return df[same_day].dropna()


def pick(site, df, key, ascending, hours, prefer_hour=None, tol=None):
    """
    Return the best window by `key` that passes the data coverage check

    :type site: str
    :param site: site name
    :type df: pandas.DataFrame
    :param df: window statistics from `rank_windows`
    :type key: str
    :param key: column to rank by
    :type ascending: bool
    :param ascending: rank smallest first if True
    :type hours: int
    :param hours: window length in hours
    :type prefer_hour: float or None
    :param prefer_hour: if given, among windows within `tol` of the best
        value, take the one whose local start hour is closest to this
    :type tol: float or None
    :param tol: tolerance on `key` used with `prefer_hour`
    :rtype: tuple or None
    :return: (window start, statistics row), or None if nothing passes
    """
    order = df.sort_values(key, ascending=ascending)
    ok = [(t, row) for t, row in order.iterrows()
          if coverage_ok(site, UTCDateTime(t.to_pydatetime()),
                         UTCDateTime(t.to_pydatetime()) + hours * 3600)]
    if not ok:
        return None
    if prefer_hour is None:
        return ok[0]
    best = ok[0][1][key]
    close = [(t, r) for t, r in ok if abs(r[key] - best) <= tol]
    tz = SITES[site]["tz"]

    def dist(item):
        """Circular distance in hours to the preferred local start hour"""
        lt = item[0].tz_convert(tz)
        d = abs(lt.hour + lt.minute / 60 - prefer_hour) % 24
        return min(d, 24 - d)

    return min(close, key=dist)


def find(site, hours=4, low_tol=0.5):
    """
    Find the high, low and change windows of one site

    :type site: str
    :param site: site name
    :type hours: int
    :param hours: window length in hours
    :type low_tol: float
    :param low_tol: m/s above the calmest window still accepted as low
        when matching the local time of the high window
    :rtype: dict
    :return: {'high'|'low'|'change': (start, statistics row)}
    """
    df = rank_windows(wind_10min(site), hours)
    out = {"high": pick(site, df, "min", False, hours)}
    lt = out["high"][0].tz_convert(SITES[site]["tz"])
    out["low"] = pick(site, df, "max", True, hours,
                      prefer_hour=lt.hour + lt.minute / 60, tol=low_tol)
    df["absdelta"] = df["delta"].abs()
    out["change"] = pick(site, df, "absdelta", False, hours)
    return out


def main():
    """Print the windows of each site as WINDOWS dict entries"""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites", nargs="+", default=list(SITES),
                        choices=list(SITES))
    parser.add_argument("--hours", type=int, default=4)
    parser.add_argument("--low_tol", type=float, default=0.5)
    args = parser.parse_args()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    for site in args.sites:
        found = find(site, args.hours, args.low_tol)
        print(f"    # {site}")
        for kind, res in found.items():
            if res is None:
                print(f"    # no {kind} window with complete data")
                continue
            t, row = res
            end = t + pd.Timedelta(hours=args.hours)
            desc = {"high": "high wind", "low": "low wind"}.get(
                kind, "rising wind" if row["delta"] > 0 else "falling wind")
            lt = t.tz_convert(SITES[site]["tz"])
            print(f'    "{site.lower()}_{kind}": ("{site}", '
                  f'{t.dayofyear}, "{t:%H:%M}", "{end:%H:%M}", "{desc}"),'
                  f"  # local {lt:%m-%d %H:%M}, mean {row['mean']:.1f} "
                  f"min {row['min']:.1f} max {row['max']:.1f} "
                  f"d {row['delta']:+.1f}")


if __name__ == "__main__":
    main()
