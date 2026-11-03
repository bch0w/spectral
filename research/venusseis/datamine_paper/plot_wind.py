#!/usr/bin/env python
"""
Plot wind roses and wind summary figures for ISHTAR anemometer stations.

Reads day-long MiniSEED files of wind speed (``LWS``, m/s) and wind
direction (``LWD``, degrees clockwise from North, direction the wind blows
FROM) from the DATA directory, merges them over the full deployment, and
averages the 1 Hz samples into fixed windows (10 min by default, the WMO
standard averaging period). Direction is averaged as a vector so that
values near 0/360 degrees do not cancel out.

Three figure types are written to FIGURES/WIND:

1. ``<STA>_windrose.png``: wind rose for each station
2. ``ALL_windrose.png``: every station's wind rose on one page with shared
   speed bins and radial scale, for direct comparison in a paper
3. ``<STA>_wind_summary.png``: wind rose plus a daily time series (median,
   interquartile range and maximum of the averaged speed) and an
   hour-of-day vs. day heatmap, which shows low- and high-wind days and
   any diurnal cycle at a glance

A CSV of daily wind statistics per station is also written so that quiet
and windy days can be cross-referenced against the seismic noise spectra.

Example::

    python plot_wind.py                       # all stations
    python plot_wind.py --stations RNDA,PAXA  # a subset
    python plot_wind.py --avg 60 --calm 0.3   # 1-min averages, 0.3 m/s calm
"""
import argparse
import os
import re
from collections import defaultdict
from glob import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from obspy import Stream, read  # noqa: E402

PATH_DATA = "./DATA"
PATH_OUT = "./FIGURES/WIND"

# Speed bin edges (m/s); last bin is open-ended
SPEED_BINS = [0.5, 2, 4, 6, 8, 10]
N_SECTORS = 16
FNAME_RE = re.compile(
    r"^(?P<net>\w+)\.(?P<sta>\w+)\.(?P<loc>\w*)\.(?P<cha>LW[DS])\."
    r"(?P<year>\d{4})\.(?P<jday>\d{3})$"
)


def find_wind_files(path_data, stations=None):
    """
    Find wind speed and direction files and group them by station.

    :type path_data: str
    :param path_data: directory containing day-long MiniSEED files
    :type stations: list of str or None
    :param stations: station codes to keep, or None to keep all
    :rtype: dict
    :return: ``{station: {"LWS": [files], "LWD": [files]}}``
    """
    files = defaultdict(lambda: {"LWS": [], "LWD": []})
    for fid in sorted(os.listdir(path_data)):
        match = FNAME_RE.match(fid)
        if not match:
            continue
        sta = match.group("sta")
        if stations and sta not in stations:
            continue
        files[sta][match.group("cha")].append(os.path.join(path_data, fid))
    return dict(files)


def read_channel(fids):
    """
    Read and merge a list of day files into a time-indexed series.

    Gaps are left out rather than filled, so they do not bias statistics.

    :type fids: list of str
    :param fids: paths to MiniSEED files for a single channel
    :rtype: pandas.Series
    :return: sample values indexed by UTC timestamp, rounded to 1 s
    """
    st = Stream()
    for fid in fids:
        st += read(fid)
    parts = []
    # Traces are not merged with ObsPy because encodings can differ
    # between days; overlaps are removed below instead
    for tr in st:
        times = (tr.stats.starttime.ns
                 + np.arange(tr.stats.npts) * int(tr.stats.delta * 1E9))
        index = pd.to_datetime(times, utc=True).round("1s")
        parts.append(pd.Series(tr.data.astype(float), index=index))
    series = pd.concat(parts).sort_index()
    return series[~series.index.duplicated(keep="first")]


def load_station(fids_lws, fids_lwd, avg_sec=600, min_frac=0.5):
    """
    Load wind speed and direction for one station and average in windows.

    Speed is the scalar mean of the samples in each window. Direction is
    the direction of the mean wind vector, so it is not biased by the
    0/360 degree wrap-around.

    :type fids_lws: list of str
    :param fids_lws: wind speed (LWS) files
    :type fids_lwd: list of str
    :param fids_lwd: wind direction (LWD) files
    :type avg_sec: int
    :param avg_sec: averaging window length in seconds; 0 or 1 keeps the
        raw 1 Hz samples
    :type min_frac: float
    :param min_frac: minimum fraction of samples a window needs to be kept
    :rtype: pandas.DataFrame
    :return: columns ``speed`` (m/s) and ``direction`` (deg), indexed by
        window start time
    """
    spd = read_channel(fids_lws)
    drc = read_channel(fids_lwd)
    df = pd.concat({"speed": spd, "direction": drc}, axis=1,
                   join="inner").dropna()
    # Remove physically impossible values
    df = df[(df.speed >= 0) & (df.speed < 100)
            & (df.direction >= 0) & (df.direction <= 360)]

    if avg_sec <= 1:
        return df

    rad = np.deg2rad(df.direction)
    # Components of the 'from' direction weighted by speed
    df = df.assign(u=df.speed * np.sin(rad), v=df.speed * np.cos(rad))
    grp = df.resample(f"{avg_sec}s")
    out = grp[["speed", "u", "v"]].mean()
    out["count"] = grp["speed"].count()
    out = out[out["count"] >= min_frac * avg_sec]
    out["direction"] = np.rad2deg(np.arctan2(out.u, out.v)) % 360
    return out[["speed", "direction"]]


def rose_histogram(df, calm=0.5, speed_bins=SPEED_BINS,
                   n_sectors=N_SECTORS):
    """
    Bin wind observations by direction sector and speed class.

    :type df: pandas.DataFrame
    :param df: wind data with ``speed`` and ``direction`` columns
    :type calm: float
    :param calm: speeds below this (m/s) are counted as calm and excluded
        from the directional bins
    :type speed_bins: list of float
    :param speed_bins: lower edges of speed classes; first must equal calm
    :type n_sectors: int
    :param n_sectors: number of direction sectors
    :rtype: tuple
    :return: (percent array of shape (n_speed, n_sectors), percent calm)
    """
    width = 360 / n_sectors
    # Center sectors on N, NNE, ... by shifting half a sector
    sector = (((df.direction + width / 2) % 360) // width).astype(int)
    edges = list(speed_bins) + [np.inf]
    moving = df.speed >= calm
    hist = np.zeros((len(speed_bins), n_sectors))
    for i in range(len(speed_bins)):
        sel = moving & (df.speed >= edges[i]) & (df.speed < edges[i + 1])
        counts = np.bincount(sector[sel], minlength=n_sectors)
        hist[i] = counts
    total = len(df)
    return 100 * hist / total, 100 * (~moving).sum() / total


def speed_labels(speed_bins=SPEED_BINS):
    """
    Build legend labels for speed classes.

    :type speed_bins: list of float
    :param speed_bins: lower edges of speed classes
    :rtype: list of str
    :return: one label per class, e.g. "2-4"
    """
    labels = [f"{a:g}–{b:g}" for a, b in zip(speed_bins[:-1],
                                                   speed_bins[1:])]
    return labels + [f"≥{speed_bins[-1]:g}"]


def plot_rose(ax, df, calm=0.5, rmax=None, title=None, legend=True,
              speed_bins=SPEED_BINS, n_sectors=N_SECTORS):
    """
    Draw a stacked wind rose on a polar axis.

    Bar length is the percentage of all observations (including calms)
    coming from that sector, stacked by speed class.

    :type ax: matplotlib.projections.polar.PolarAxes
    :param ax: polar axis to draw on
    :type df: pandas.DataFrame
    :param df: wind data with ``speed`` and ``direction`` columns
    :type calm: float
    :param calm: calm threshold in m/s
    :type rmax: float or None
    :param rmax: radial axis limit in percent; None scales to the data
    :type title: str or None
    :param title: axis title
    :type legend: bool
    :param legend: draw a speed legend next to the axis
    :type speed_bins: list of float
    :param speed_bins: lower edges of speed classes
    :type n_sectors: int
    :param n_sectors: number of direction sectors
    :rtype: float
    :return: largest sector total in percent
    """
    hist, pct_calm = rose_histogram(df, calm, speed_bins, n_sectors)
    cmap = plt.get_cmap("viridis", len(speed_bins))
    theta = np.deg2rad(np.arange(n_sectors) * 360 / n_sectors)
    width = np.deg2rad(360 / n_sectors) * 0.9
    bottom = np.zeros(n_sectors)
    for i, label in enumerate(speed_labels(speed_bins)):
        ax.bar(theta, hist[i], width=width, bottom=bottom,
               color=cmap(i), edgecolor="white", linewidth=0.4,
               label=label)
        bottom += hist[i]

    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_xticks(np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315]))
    ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
    peak = bottom.max()
    rlim = rmax or peak * 1.1
    ax.set_ylim(0, rlim)
    # Keep radial labels sparse (about 4 rings) so they stay legible
    step = max(1, np.ceil(rlim / 4 / 5) * 5) if rlim > 8 else 2
    ax.set_yticks(np.arange(step, rlim, step))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda r, _: f"{r:g}%"))
    ax.set_rlabel_position(112.5)
    ax.tick_params(axis="y", labelsize=7)
    ax.grid(alpha=0.4)
    ax.text(0.5, -0.18, f"calm (<{calm:g} m/s): {pct_calm:.1f}%",
            transform=ax.transAxes, ha="center", fontsize=8)
    if title:
        ax.set_title(title, pad=14)
    if legend:
        ax.legend(title="Speed (m/s)", loc="center left",
                  bbox_to_anchor=(1.1, 0.5), fontsize=8, title_fontsize=8,
                  frameon=False)
    return peak


def date_span(df):
    """
    Format the date range of a dataset for figure titles.

    :type df: pandas.DataFrame
    :param df: time-indexed data
    :rtype: str
    :return: e.g. "2025-03-26 to 2025-05-07 (42 d)"
    """
    t0, t1 = df.index[0], df.index[-1]
    ndays = (t1 - t0).total_seconds() / 86400
    return f"{t0:%Y-%m-%d} to {t1:%Y-%m-%d} ({ndays:.0f} d)"


def daily_stats(df, calm=0.5):
    """
    Compute daily wind speed statistics from averaged wind data.

    :type df: pandas.DataFrame
    :param df: averaged wind data with ``speed`` and ``direction`` columns
    :type calm: float
    :param calm: calm threshold in m/s
    :rtype: pandas.DataFrame
    :return: one row per UTC day with median, percentiles, max, percent
        calm and data coverage
    """
    grp = df.speed.resample("1D")
    stats = pd.DataFrame({
        "median": grp.median(),
        "p25": grp.quantile(0.25),
        "p75": grp.quantile(0.75),
        "p95": grp.quantile(0.95),
        "max": grp.max(),
        "pct_calm": grp.apply(lambda s: 100 * (s < calm).mean()
                              if len(s) else np.nan),
        "n": grp.count(),
    })
    return stats.dropna(subset=["median"])


def plot_summary(df, sta, fid_out, calm=0.5, avg_sec=600):
    """
    Plot a three-panel wind summary for one station.

    Panels: wind rose; daily median speed with interquartile band and
    daily maximum; hour-of-day vs. date heatmap of mean speed.

    :type df: pandas.DataFrame
    :param df: averaged wind data with ``speed`` and ``direction`` columns
    :type sta: str
    :param sta: station code for the title
    :type fid_out: str
    :param fid_out: output figure path
    :type calm: float
    :param calm: calm threshold in m/s
    :type avg_sec: int
    :param avg_sec: averaging window used, for labeling
    """
    stats = daily_stats(df, calm)
    fig = plt.figure(figsize=(13, 7.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 2.1],
                          height_ratios=[1, 1.1], wspace=0.25, hspace=0.35)
    ax_rose = fig.add_subplot(gs[:, 0], projection="polar")
    ax_ts = fig.add_subplot(gs[0, 1])
    ax_hm = fig.add_subplot(gs[1, 1], sharex=ax_ts)

    plot_rose(ax_rose, df, calm=calm, legend=False)
    ax_rose.legend(title="Speed (m/s)", loc="upper center",
                   bbox_to_anchor=(0.5, -0.18), ncol=3, fontsize=8,
                   title_fontsize=8, frameon=False)

    # Daily time series
    days = stats.index + pd.Timedelta(hours=12)
    # Every averaged sample at its true time, faint, to show the spread
    ax_ts.scatter(df.index + pd.Timedelta(seconds=avg_sec / 2), df.speed,
                  s=2, color="0.35", alpha=0.15, lw=0, zorder=0,
                  rasterized=True, label=f"{avg_sec // 60:g}-min mean")
    ax_ts.fill_between(days, stats.p25, stats.p75, color="C0", alpha=0.3,
                       lw=0, label="interquartile range")
    ax_ts.plot(days, stats["median"], "o-", color="C0", ms=3,
               label="median")
    ax_ts.plot(days, stats["max"], "v", color="C3", ms=4,
               label=f"max {avg_sec // 60:g}-min mean")
    q_lo, q_hi = stats["median"].quantile([0.25, 0.75])
    ax_ts.axhspan(0, q_lo, color="0.5", alpha=0.08, lw=0)
    ax_ts.axhline(q_hi, color="0.4", lw=0.8, ls="--")
    ax_ts.axhline(q_lo, color="0.4", lw=0.8, ls=":")
    ax_ts.text(1.005, q_hi, "windy-day\nthreshold", va="center",
               fontsize=7, transform=ax_ts.get_yaxis_transform())
    ax_ts.text(1.005, q_lo, "quiet-day\nthreshold", va="center",
               fontsize=7, transform=ax_ts.get_yaxis_transform())
    ax_ts.set_ylabel("Daily wind speed (m/s)")
    ax_ts.set_ylim(bottom=0)
    leg = ax_ts.legend(loc="lower left", bbox_to_anchor=(0, 1.0),
                       fontsize=8, ncol=4, frameon=False,
                       borderaxespad=0.2)
    # Legend marker for the faint scatter is otherwise nearly invisible
    for handle in leg.legend_handles:
        if hasattr(handle, "set_sizes"):
            handle.set_sizes([12])
            handle.set_alpha(0.6)
    ax_ts.grid(alpha=0.3)
    plt.setp(ax_ts.get_xticklabels(), visible=False)

    # Hour-of-day vs date heatmap (hourly mean speed)
    hourly = df.speed.resample("1h").mean()
    grid = hourly.to_frame("speed").assign(
        date=hourly.index.floor("D"), hour=hourly.index.hour
    ).pivot_table(index="hour", columns="date", values="speed")
    grid = grid.reindex(columns=pd.date_range(grid.columns[0],
                                              grid.columns[-1], freq="D"))
    x = mdates.date2num(list(grid.columns) + [grid.columns[-1]
                                               + pd.Timedelta(days=1)])
    mesh = ax_hm.pcolormesh(x, np.arange(25), grid.values, cmap="viridis",
                            shading="flat", vmin=0)
    cax = ax_hm.inset_axes([1.01, 0, 0.015, 1])
    fig.colorbar(mesh, cax=cax, label="Hourly mean speed (m/s)")
    ax_hm.set_ylabel("Hour of day (UTC)")
    ax_hm.set_yticks([0, 6, 12, 18, 24])
    ax_hm.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax_hm.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))

    fig.suptitle(f"{sta} wind, {date_span(df)}", fontsize=13)
    fig.savefig(fid_out, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_all_roses(data, fid_out, calm=0.5):
    """
    Plot every station's wind rose on one figure with a shared scale.

    :type data: dict
    :param data: ``{station: DataFrame}`` of averaged wind data
    :type fid_out: str
    :param fid_out: output figure path
    :type calm: float
    :param calm: calm threshold in m/s
    """
    # Order panels by deployment start time, then station code for ties
    stas = sorted(data, key=lambda s: (data[s].index[0], s))
    peaks = [rose_histogram(data[s], calm)[0].sum(axis=0).max()
             for s in stas]
    rmax = np.ceil(max(peaks) * 1.05)
    ncol = min(4, len(stas))
    nrow = int(np.ceil(len(stas) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 4.4 * nrow),
                            subplot_kw={"projection": "polar"},
                            squeeze=False,
                            gridspec_kw={"hspace": 0.55, "wspace": 0.35})
    for ax, sta in zip(axs.flat, stas):
        plot_rose(ax, data[sta], calm=calm, rmax=rmax, legend=False,
                  title=f"{sta}\n{date_span(data[sta])}")
        ax.title.set_fontsize(9)
    for ax in axs.flat[len(stas):]:
        ax.set_visible(False)
    handles, labels = axs.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, title="Speed (m/s)", loc="lower center",
               ncol=len(labels), frameon=False,
               bbox_to_anchor=(0.5, -0.02 - 0.06 / nrow))
    fig.savefig(fid_out, dpi=200, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    """
    Parse command line arguments.

    :rtype: argparse.Namespace
    :return: parsed arguments
    """
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data", default=PATH_DATA,
                        help="directory of MiniSEED day files")
    parser.add_argument("--out", default=PATH_OUT,
                        help="output figure directory")
    parser.add_argument("--stations", default=None,
                        help="comma-separated station codes, e.g. RNDA,PAXA")
    parser.add_argument("--avg", type=int, default=600,
                        help="averaging window in seconds (default 600; "
                             "1 uses raw 1 Hz samples)")
    parser.add_argument("--calm", type=float, default=SPEED_BINS[0],
                        help="calm threshold in m/s (default 0.5)")
    return parser.parse_args()


def main():
    """
    Read wind data for every station and write all wind figures.
    """
    args = parse_args()
    stations = args.stations.split(",") if args.stations else None
    os.makedirs(args.out, exist_ok=True)
    SPEED_BINS[0] = args.calm

    files = find_wind_files(args.data, stations)
    data = {}
    for sta, fids in sorted(files.items()):
        if not fids["LWS"] or not fids["LWD"]:
            print(f"{sta}: missing LWS or LWD files, skipping")
            continue
        print(f"{sta}: reading {len(fids['LWS'])} LWS and "
              f"{len(fids['LWD'])} LWD files")
        df = load_station(fids["LWS"], fids["LWD"], avg_sec=args.avg)
        if df.empty:
            print(f"{sta}: no overlapping speed/direction data, skipping")
            continue
        data[sta] = df
        print(f"{sta}: {len(df)} samples, {date_span(df)}, "
              f"median {df.speed.median():.2f} m/s, "
              f"max {df.speed.max():.2f} m/s")

        fig, ax = plt.subplots(figsize=(6.5, 5.5),
                               subplot_kw={"projection": "polar"})
        plot_rose(ax, df, calm=args.calm,
                  title=f"{sta} wind rose\n{date_span(df)}")
        fig.savefig(os.path.join(args.out, f"{sta}_windrose.png"), dpi=200,
                    bbox_inches="tight")
        plt.close(fig)

        plot_summary(df, sta, os.path.join(args.out,
                                           f"{sta}_wind_summary.png"),
                     calm=args.calm, avg_sec=args.avg)
        daily_stats(df, args.calm).round(3).to_csv(
            os.path.join(args.out, f"{sta}_daily_wind.csv"),
            index_label="date"
        )

    if len(data) > 1:
        plot_all_roses(data, os.path.join(args.out, "ALL_windrose.png"),
                       calm=args.calm)
    print(f"figures written to {args.out}")


if __name__ == "__main__":
    main()
