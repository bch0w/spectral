#!/usr/bin/env python
"""
Long-term spectrograms to check how common a narrowband signal is.

For each channel, every day file is cut into segments (10 min by default),
the power spectral density of each segment is estimated with Welch's method
and converted to velocity with the scalar instrument sensitivity (fine
above the sensor corner frequency), and the results are cached in an .npz
file. Two figures are made from the cache:

- ``<name>_time.png``: frequency vs. time over the whole deployment, one
  panel per channel, with the anemometer wind speed on top
- ``<name>_tod.png``: median PSD per local hour of day (frequency vs. hour),
  one panel per channel, to show diurnal (e.g. temperature driven) patterns

Example::

    python longterm_psd.py XX.DESDC..HHN XX.DESDU..HHN XX.DESDW..HHN \\
        HV.DESD..EHN --wind XX.DESDA..LWS --tz Pacific/Honolulu \\
        --name desd_N --path_out CHECKS/desd_35hz
"""
import argparse
import glob
import os
from datetime import datetime, timedelta
from zoneinfo import ZoneInfo

import matplotlib
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from obspy import read, read_inventory, UTCDateTime
from scipy.signal import welch

PATH_DATA = "DATA"
PATH_XML = "XML"


def sensitivity(code, time, inv):
    """
    Return the scalar sensitivity (counts per m/s) of a channel

    :type code: str
    :param code: trace ID, e.g. 'XX.DESDU..HHN'
    :type time: obspy.UTCDateTime
    :param time: time at which to look up the response
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory holding the channel
    :rtype: float
    :return: instrument sensitivity
    """
    resp = inv.get_response(code, time)
    return resp.instrument_sensitivity.value


def channel_psds(code, inv, seg=600., nperseg=1024, min_frac=0.9):
    """
    Compute the PSD of every segment of every day file of one channel

    :type code: str
    :param code: trace ID, e.g. 'XX.DESDU..HHN'
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory holding the channel
    :type seg: float
    :param seg: segment length in seconds
    :type nperseg: int
    :param nperseg: Welch sub-window length in samples
    :type min_frac: float
    :param min_frac: minimum data fraction a segment needs to be used
    :rtype: tuple of np.array
    :return: segment start times (POSIX seconds), frequencies (Hz) and
        PSDs in dB re (m/s)^2/Hz with shape (time, frequency)
    """
    times, psds, freq = [], [], None
    cha = code.split(".")[-1]
    for fid in sorted(glob.glob(os.path.join(PATH_DATA, f"{code}.*"))):
        st = read(fid).select(channel=cha)
        if not st:
            continue
        st.merge(method=1, fill_value=None)
        tr = st[0]
        sens = sensitivity(code, tr.stats.starttime, inv)
        fs = tr.stats.sampling_rate
        day = UTCDateTime(tr.stats.starttime.date)
        nseg = int(86400 // seg)
        for i in range(nseg):
            t0 = day + i * seg
            sl = tr.slice(t0, t0 + seg - tr.stats.delta)
            data = sl.data
            if np.ma.isMaskedArray(data):
                data = data.compressed()
            if len(data) < min_frac * seg * fs:
                continue
            f, p = welch(data.astype(float) / sens, fs=fs,
                         nperseg=nperseg, detrend="linear")
            freq = f
            times.append(t0.timestamp)
            psds.append(10 * np.log10(np.maximum(p, 1E-30)))
        print(f"   {os.path.basename(fid)}: {len(times)} segments so far")
    return np.array(times), freq, np.array(psds)


def load_or_compute(code, inv, path_cache, seg, nperseg):
    """
    Return the cached PSDs of a channel, computing them if needed

    :type code: str
    :param code: trace ID
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory holding the channel
    :type path_cache: str
    :param path_cache: directory of the .npz cache files
    :type seg: float
    :param seg: segment length in seconds
    :type nperseg: int
    :param nperseg: Welch sub-window length in samples
    :rtype: tuple of np.array
    :return: times, frequencies and PSDs as from `channel_psds`
    """
    fid = os.path.join(path_cache, f"{code}_{seg:g}s_{nperseg}.npz")
    if os.path.exists(fid):
        d = np.load(fid)
        return d["times"], d["freq"], d["psds"]
    print(f"{code}: computing PSDs")
    times, freq, psds = channel_psds(code, inv, seg, nperseg)
    os.makedirs(path_cache, exist_ok=True)
    np.savez_compressed(fid, times=times, freq=freq, psds=psds)
    return times, freq, psds


def gridded(times, psds, seg):
    """
    Put PSDs on a regular time grid, with NaN where data are missing, so
    gaps show as blank columns

    :type times: np.array
    :param times: segment start times (POSIX seconds)
    :type psds: np.array
    :param psds: PSDs, shape (time, frequency)
    :type seg: float
    :param seg: segment length in seconds
    :rtype: tuple of np.array
    :return: regular start times and PSDs with NaN rows for gaps
    """
    grid = np.arange(times.min(), times.max() + seg, seg)
    out = np.full((len(grid), psds.shape[1]), np.nan)
    out[np.round((times - grid[0]) / seg).astype(int)] = psds
    return grid, out


def wind_series(code, seg):
    """
    Read every wind speed day file and average it in segments

    :type code: str
    :param code: wind speed trace ID, e.g. 'XX.DESDA..LWS'
    :type seg: float
    :param seg: averaging length in seconds
    :rtype: tuple of np.array
    :return: segment start times (POSIX seconds) and mean speed (m/s);
        segments holding fault values (<0 or >=100 m/s) are dropped
    """
    times, speeds = [], []
    for fid in sorted(glob.glob(os.path.join(PATH_DATA, f"{code}.*"))):
        for tr in read(fid):
            t = tr.stats.starttime.timestamp + tr.times()
            times.append(t)
            speeds.append(tr.data.astype(float))
    t, s = np.concatenate(times), np.concatenate(speeds)
    bins = np.floor(t / seg) * seg
    ubins, idx = np.unique(bins, return_inverse=True)
    mean = np.bincount(idx, weights=s) / np.bincount(idx)
    bad = np.bincount(idx, weights=((s < 0) | (s >= 100)).astype(float))
    return ubins[bad == 0], mean[bad == 0]


def plot_time(results, wind, fmin, fmax, vmin, vmax, tz, fid_out):
    """
    Plot frequency vs. time for each channel over the whole record

    :type results: dict
    :param results: {code: (times, freq, psds)}
    :type wind: tuple or None
    :param wind: (times, speed) of the wind, or None
    :type fmin: float
    :param fmin: lowest frequency shown
    :type fmax: float
    :param fmax: highest frequency shown
    :type vmin: float
    :param vmin: lower color limit in dB
    :type vmax: float
    :param vmax: upper color limit in dB
    :type tz: str
    :param tz: time zone for the local-midnight markers
    :type fid_out: str
    :param fid_out: output figure path
    """
    n = len(results) + (1 if wind else 0)
    ratios = ([0.5] if wind else []) + [1] * len(results)
    f, axs = plt.subplots(n, 1, figsize=(16, 2.6 * n), sharex=True,
                          gridspec_kw={"height_ratios": ratios})
    axs = np.atleast_1d(axs)
    i = 0
    if wind:
        wt, ws = wind
        # Matplotlib date numbers are days since 1970-01-01
        axs[0].plot(wt / 86400., ws, c="C3", lw=0.5)
        axs[0].set_ylabel("Wind (m/s)\n10-min mean")
        i = 1
    t_all = []
    for ax, (code, (times, freq, psds)) in zip(axs[i:], results.items()):
        seg = np.median(np.diff(times))
        grid, data = gridded(times, psds, seg)
        keep = (freq >= fmin) & (freq <= fmax)
        x = grid / 86400.  # Matplotlib date numbers (days since epoch)
        im = ax.pcolormesh(np.append(x, x[-1] + seg / 86400.),
                           freq[keep], data[:, keep].T[:-1],
                           vmin=vmin, vmax=vmax, cmap="viridis",
                           shading="flat", rasterized=True)
        # The temporary XX network code is left out of figure labels
        ax.set_ylabel(f"{code.removeprefix('XX.')}\nFreq. (Hz)")
        t_all.append(x)
        f.colorbar(im, ax=ax, pad=0.01, label="dB re (m/s)$^2$/Hz")
    if wind:
        # Same width as the panels below, but no visible colorbar
        f.colorbar(im, ax=axs[0], pad=0.01).ax.set_visible(False)
    # Local midnight markers
    x = np.concatenate(t_all)
    day = UTCDateTime(x.min() * 86400).date
    zone = ZoneInfo(tz)
    d = datetime(day.year, day.month, day.day, tzinfo=zone)
    while d.timestamp() < x.max() * 86400:
        for ax in axs:
            ax.axvline(d.timestamp() / 86400., c="w" if ax is not axs[0]
                       or not wind else "0.7", lw=0.4, alpha=0.6)
        d += timedelta(days=1)
    axs[-1].xaxis.set_major_locator(mdates.DayLocator(interval=2))
    axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%j"))
    axs[-1].set_xlabel(f"Date (UTC, month-day / julian day); lines at "
                       f"local midnight ({tz})")
    axs[-1].set_xlim(x.min(), x.max())
    f.tight_layout()
    f.savefig(fid_out, dpi=150)
    plt.close(f)


def plot_tod(results, fmin, fmax, vmin, vmax, tz, fid_out):
    """
    Plot the median PSD per local hour of day for each channel

    :type results: dict
    :param results: {code: (times, freq, psds)}
    :type fmin: float
    :param fmin: lowest frequency shown
    :type fmax: float
    :param fmax: highest frequency shown
    :type vmin: float
    :param vmin: lower color limit in dB
    :type vmax: float
    :param vmax: upper color limit in dB
    :type tz: str
    :param tz: time zone for the hour of day
    :type fid_out: str
    :param fid_out: output figure path
    """
    zone = ZoneInfo(tz)
    n = len(results)
    f, axs = plt.subplots(1, n, figsize=(4.2 * n, 4), sharey=True)
    axs = np.atleast_1d(axs)
    for ax, (code, (times, freq, psds)) in zip(axs, results.items()):
        hours = np.array([
            UTCDateTime(t).datetime.replace(
                tzinfo=ZoneInfo("UTC")).astimezone(zone).hour
            for t in times])
        keep = (freq >= fmin) & (freq <= fmax)
        med = np.full((24, keep.sum()), np.nan)
        for h in range(24):
            if (hours == h).any():
                med[h] = np.median(psds[hours == h][:, keep], axis=0)
        im = ax.pcolormesh(np.arange(25), freq[keep], med.T[:-1],
                           vmin=vmin, vmax=vmax, cmap="viridis",
                           shading="flat")
        ax.set_title(f"{code.removeprefix('XX.')}\n"
                     f"({len(times)} segments, "
                     f"{len(set(np.floor(times / 86400)))} days)",
                     fontsize=10)
        ax.set_xlabel(f"Local hour ({tz})")
        ax.set_xticks(range(0, 25, 6))
    axs[0].set_ylabel("Freq. (Hz)")
    f.colorbar(im, ax=axs, pad=0.01, label="median dB re (m/s)$^2$/Hz")
    f.savefig(fid_out, dpi=150, bbox_inches="tight")
    plt.close(f)


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("codes", nargs="+",
                        help="trace IDs, e.g. XX.DESDU..HHN")
    parser.add_argument("--wind", default=None,
                        help="wind speed trace ID to plot on top, e.g. "
                             "XX.DESDA..LWS")
    parser.add_argument("--tz", default="UTC",
                        help="local time zone, e.g. Pacific/Honolulu")
    parser.add_argument("--seg", type=float, default=600.,
                        help="segment length in seconds, default 600")
    parser.add_argument("--nperseg", type=int, default=1024,
                        help="Welch sub-window in samples, default 1024")
    parser.add_argument("--fmin", type=float, default=5.)
    parser.add_argument("--fmax", type=float, default=50.)
    parser.add_argument("--vmin", type=float, default=-190.)
    parser.add_argument("--vmax", type=float, default=-140.)
    parser.add_argument("--name", default="longterm",
                        help="figure name prefix")
    parser.add_argument("--path_out", default="CHECKS")
    return parser.parse_args()


def main():
    """Compute or load PSDs for each channel and make both figures"""
    args = parse_args()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    matplotlib.use("Agg")
    inv = read_inventory(os.path.join(PATH_XML, "ISHTAR.xml"))
    for code in args.codes:
        if not code.startswith("XX."):
            netsta = ".".join(code.split(".")[:2])
            inv += read_inventory(os.path.join(PATH_XML, f"{netsta}.xml"))
    path_cache = os.path.join(args.path_out, "psd_cache")
    results = {code: load_or_compute(code, inv, path_cache, args.seg,
                                     args.nperseg)
               for code in args.codes}
    wind = wind_series(args.wind, args.seg) if args.wind else None
    os.makedirs(args.path_out, exist_ok=True)
    fid = os.path.join(args.path_out, f"{args.name}_time.png")
    plot_time(results, wind, args.fmin, args.fmax, args.vmin, args.vmax,
              args.tz, fid)
    print(f"-> {fid}")
    fid = os.path.join(args.path_out, f"{args.name}_tod.png")
    plot_tod(results, args.fmin, args.fmax, args.vmin, args.vmax, args.tz,
             fid)
    print(f"-> {fid}")


if __name__ == "__main__":
    main()
