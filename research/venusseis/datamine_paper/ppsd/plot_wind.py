#!/usr/bin/env python
"""
Plot wind roses and wind speed figures for ISHTAR anemometer stations.

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
3. ``ALL_wind_speed.png``: average wind speed over time for every station
   with one column per deployment period. Each averaged sample is colored
   by direction, with the daily median, interquartile range and a
   downwind arrow at each day's maximum, which shows low- and high-wind
   days at a glance

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
# Cyclic colormap for wind direction: N red, E gold, S teal, W blue, and
# back to red, so 359 and 1 degrees share a color. Colors are kept at mid
# lightness so every direction stays visible on a white background
DIR_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list(
    "wind_direction",
    ["#c0392b", "#e09f1f", "#2a9d8f", "#3b5ba5", "#c0392b"], N=256
)
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
    # Fault codes (e.g. speed 1000, direction 999) come in bursts where the
    # neighboring 'valid' samples are also corrupt (20-100 m/s), so every
    # window containing a fault is dropped, not just the fault samples
    fault = ((df.speed < 0) | (df.speed >= 100)
             | (df.direction < 0) | (df.direction > 360))
    block = f"{max(avg_sec, 60)}s"
    bad_blocks = df.index[fault].floor(block).unique()
    in_bad = df.index.floor(block).isin(bad_blocks)
    if fault.any():
        print(f"  dropping {len(bad_blocks)} {block} window(s) with "
              f"{fault.sum()} fault samples")
    df = df[~in_bad]

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
              speed_bins=SPEED_BINS, n_sectors=N_SECTORS, fontsize=10,
              calm_y=-0.18):
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
    :type fontsize: float
    :param fontsize: font size of the compass labels, calm note and legend;
        radial labels are 3 pt smaller
    :type calm_y: float
    :param calm_y: vertical position of the calm note in axis coordinates,
        or None to leave the note out
    :rtype: float
    :return: largest sector total in percent
    """
    hist, pct_calm = rose_histogram(df, calm, speed_bins, n_sectors)
    cmap = plt.get_cmap("plasma", len(speed_bins))
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
    ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"],
                       fontsize=fontsize)
    peak = bottom.max()
    rlim = rmax or peak * 1.1
    ax.set_ylim(0, rlim)
    # Keep radial labels sparse (about 4 rings) so they stay legible
    step = max(1, np.ceil(rlim / 4 / 5) * 5) if rlim > 8 else 2
    ax.set_yticks(np.arange(step, rlim, step))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda r, _: f"{r:g}%"))
    ax.set_rlabel_position(112.5)
    ax.tick_params(axis="y", labelsize=fontsize - 3)
    ax.grid(alpha=0.4)
    if calm_y is not None:
        ax.text(0.5, calm_y, f"calm (<{calm:g} m/s): {pct_calm:.1f}%",
                transform=ax.transAxes, ha="center", fontsize=fontsize)
    if title:
        ax.set_title(title, pad=14)
    if legend:
        ax.legend(title="Speed (m/s)", loc="center left",
                  bbox_to_anchor=(1.1, 0.5), fontsize=fontsize,
                  title_fontsize=fontsize,
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


def date_span_short(df):
    """
    Format the date range of a dataset compactly for small panel titles.

    :type df: pandas.DataFrame
    :param df: time-indexed data
    :rtype: str
    :return: e.g. "2025-03-26 – 05-07 (42 d)"; the end year is only
        repeated if it differs from the start year
    """
    t0, t1 = df.index[0], df.index[-1]
    ndays = (t1 - t0).total_seconds() / 86400
    end = f"{t1:%m-%d}" if t1.year == t0.year else f"{t1:%Y-%m-%d}"
    return f"{t0:%Y-%m-%d} – {end} ({ndays:.0f} d)"


def daily_stats(df, calm=0.5):
    """
    Compute daily wind speed statistics from averaged wind data.

    :type df: pandas.DataFrame
    :param df: averaged wind data with ``speed`` and ``direction`` columns
    :type calm: float
    :param calm: calm threshold in m/s
    :rtype: pandas.DataFrame
    :return: one row per UTC day with median, percentiles, max, direction
        (from) at the time of the max, percent calm and data coverage
    """
    grp = df.speed.resample("1D")
    # Time of each day's strongest averaged wind, for its direction
    t_max = grp.apply(lambda s: s.idxmax() if len(s) else pd.NaT)
    dir_at_max = pd.Series(
        [df.direction.get(t, np.nan) if pd.notna(t) else np.nan
         for t in t_max], index=t_max.index
    )
    stats = pd.DataFrame({
        "median": grp.median(),
        "p25": grp.quantile(0.25),
        "p75": grp.quantile(0.75),
        "p95": grp.quantile(0.95),
        "max": grp.max(),
        "dir_at_max": dir_at_max,
        "pct_calm": grp.apply(lambda s: 100 * (s < calm).mean()
                              if len(s) else np.nan),
        "n": grp.count(),
    })
    return stats.dropna(subset=["median"])


def plot_direction_wheel(ax, fontsize=9):
    """
    Draw a color ring on a polar axis as the legend for DIR_CMAP.

    :type ax: matplotlib.projections.polar.PolarAxes
    :param ax: polar axis to draw the ring on
    :type fontsize: float
    :param fontsize: font size of the compass labels and title
    """
    theta = np.linspace(0, 2 * np.pi, 361)
    radius = np.array([0.55, 1.0])
    color = np.tile(np.rad2deg(theta[:-1]), (1, 1))
    ax.pcolormesh(theta, radius, color, cmap=DIR_CMAP, vmin=0, vmax=360,
                  shading="flat", rasterized=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks(np.deg2rad([0, 90, 180, 270]))
    ax.set_xticklabels(["N", "E", "S", "W"], fontsize=fontsize)
    ax.tick_params(axis="x", pad=-1)
    ax.grid(False)
    ax.spines["polar"].set_visible(False)
    ax.set_title("Direction (from)", fontsize=fontsize, pad=4)


def plot_speed_panel(ax, df, calm=0.5, avg_sec=600, fontsize=10,
                     ylabel=True):
    """
    Plot averaged wind speed over time for one station on one axis.

    Every averaged sample is drawn at its true time and colored by
    direction (calms in gray). The daily median and interquartile band are
    overlaid, and each day's maximum is an arrow pointing downwind.

    :type ax: matplotlib.axes.Axes
    :param ax: axis to draw on
    :type df: pandas.DataFrame
    :param df: averaged wind data with ``speed`` and ``direction`` columns
    :type calm: float
    :param calm: calm threshold in m/s
    :type avg_sec: int
    :param avg_sec: averaging window used, for labeling
    :type fontsize: float
    :param fontsize: font size of tick and axis labels
    :type ylabel: bool
    :param ylabel: label the y-axis; turn off when the figure carries a
        single shared y-label
    :rtype: float
    :return: upper y-limit used, with headroom for the arrows
    """
    stats = daily_stats(df, calm)
    days = stats.index + pd.Timedelta(hours=12)
    times = df.index + pd.Timedelta(seconds=avg_sec / 2)
    is_calm = (df.speed < calm).values
    ax.scatter(times[is_calm], df.speed[is_calm], s=2.5, color="0.6",
               alpha=0.4, lw=0, zorder=0, rasterized=True)
    ax.scatter(times[~is_calm], df.speed[~is_calm],
               c=df.direction[~is_calm], cmap=DIR_CMAP, vmin=0, vmax=360,
               s=2.5, alpha=0.5, lw=0, zorder=0, rasterized=True,
               label=f"{avg_sec // 60:g}-min mean")
    # Summary overlays are neutral so they do not compete with the
    # direction colors
    ax.fill_between(days, stats.p25, stats.p75, color="0.2", alpha=0.15,
                    lw=0)
    # Faint dashed edges so the band reads at print size
    for edge in (stats.p25, stats.p75):
        ax.plot(days, edge, ls="--", lw=0.8, color="0.3", alpha=0.6,
                zorder=1)
    # Legend proxy showing both the fill and its dashed edge
    ax.add_patch(matplotlib.patches.Rectangle(
        (0, 0), 0, 0, facecolor=(0.2, 0.2, 0.2, 0.15),
        edgecolor=(0.3, 0.3, 0.3, 0.6), ls="--", lw=0.8, label="daily IQR"
    ))
    ax.plot(days, stats["median"], "o-", color="k", ms=2.5, lw=1,
            label="daily median")
    # Daily max as an arrow pointing downwind (direction the wind blows
    # toward) in map view: north up, east right. 'uv' angles are in screen
    # space, so orientation does not depend on the axis scaling
    to_rad = np.deg2rad((stats.dir_at_max + 180) % 360)
    ax.quiver(mdates.date2num(days), stats["max"], np.sin(to_rad),
              np.cos(to_rad), angles="uv", pivot="middle",
              scale_units="inches", scale=6, units="inches", width=0.01,
              headwidth=4, headlength=4.5, headaxislength=4, color="k",
              zorder=3)
    ax.plot([], [], ls="none", marker=r"$\rightarrow$", ms=fontsize + 3,
            color="k", label="daily max\n(downwind)")
    if ylabel:
        ax.set_ylabel("Average wind speed (m/s)", fontsize=fontsize)
    # Headroom so arrows at the top are not clipped
    ytop = stats["max"].max() * 1.12
    ax.set_ylim(0, ytop)
    ax.grid(alpha=0.3)
    ax.tick_params(labelsize=fontsize)
    # Round tick steps keep labels short (no 2.5/12.5), which matters when
    # panels sit close together
    ax.yaxis.set_major_locator(
        matplotlib.ticker.MaxNLocator(nbins=5, steps=[1, 2, 5, 10])
    )
    ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=3,
                                                      maxticks=5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    return ytop


def plot_speed_overview(data, fid_out, calm=0.5, avg_sec=600,
                        splits=("2025-07-01", "2025-08-15"),
                        own_ylim=("K24KA",)):
    """
    Plot every station's average wind speed panel, one column per campaign.

    Deployments are grouped into columns by start time using the
    ``splits`` dates, so stations in a column cover similar date ranges;
    each column is ordered by start time. Panels are labeled
    alphabetically down each column in turn. A single direction color
    wheel and the shared legend fill the empty cells of the column with
    the most free space.

    :type data: dict
    :param data: ``{station: DataFrame}`` of averaged wind data
    :type fid_out: str
    :param fid_out: output figure path
    :type calm: float
    :param calm: calm threshold in m/s
    :type avg_sec: int
    :param avg_sec: averaging window used, for labeling
    :type splits: tuple of str
    :param splits: UTC dates separating the columns, in increasing order;
        N dates give N + 1 columns
    :type own_ylim: tuple of str
    :param own_ylim: stations that keep their own y-axis range; all other
        panels in a column share one range so speeds compare directly
    """
    edges = ([pd.Timestamp.min.tz_localize("UTC")]
             + [pd.Timestamp(s, tz="UTC") for s in splits]
             + [pd.Timestamp.max.tz_localize("UTC")])
    order = sorted(data, key=lambda s: (data[s].index[0], s))
    cols = [[s for s in order if lo <= data[s].index[0] < hi]
            for lo, hi in zip(edges[:-1], edges[1:])]
    cols = [c for c in cols if c]
    if not cols:
        return
    ncol = len(cols)
    nrow = max(len(c) for c in cols)

    # Sized close to a full-page journal width (~7 in) so that, once
    # scaled down, the smallest text prints at roughly 6-7 pt
    fs = 10
    fig = plt.figure(figsize=(3.9 * ncol, 2.5 * nrow + 0.4))
    # Tight spacing: x tick labels are shown only on the bottom panel of
    # each column, so rows need room for the panel titles alone
    gs = fig.add_gridspec(nrow, ncol, hspace=0.2, wspace=0.13)
    letters = iter("abcdefghijklmnopqrstuvwxyz")
    first_ax = None
    for icol, stas in enumerate(cols):
        # Every panel in a column spans the same whole UTC days, from the
        # earliest start to the latest end in that column, so dates line up
        # vertically even where a station has no data
        t0 = min(data[s].index[0] for s in stas).floor("D")
        t1 = (max(data[s].index[-1] for s in stas)
              + pd.Timedelta(seconds=avg_sec)).ceil("D")
        col_ax = None
        shared = {}
        for irow, sta in enumerate(stas):
            ax = fig.add_subplot(gs[irow, icol], sharex=col_ax)
            col_ax = col_ax or ax
            ytop = plot_speed_panel(ax, data[sta], calm=calm,
                                    avg_sec=avg_sec, fontsize=fs,
                                    ylabel=False)
            if sta not in own_ylim:
                shared[ax] = ytop
            ax.set_title(f"({next(letters)}) {sta}", loc="left",
                         fontsize=fs + 2, fontweight="bold")
            # Compact span; full dates are on the shared time axis
            t_a, t_b = data[sta].index[0], data[sta].index[-1]
            ndays = (t_b - t_a).total_seconds() / 86400
            ax.set_title(f"{t_a:%b %d}–{t_b:%b %d} ({ndays:.0f} d)",
                         loc="right", fontsize=fs - 1)
            first_ax = first_ax or ax
        # Half-day pad keeps arrows on the first and last days unclipped
        pad = pd.Timedelta(hours=12)
        col_ax.set_xlim(t0 - pad, t1 + pad)
        # Date labels only under the bottom panel of the column
        for ax in fig.axes[-len(stas):-1]:
            ax.tick_params(labelbottom=False)
        # Evenly spaced date ticks from the column start; the automatic
        # locator gives only two ticks for month-long spans at this width
        ndays = (t1 - t0).days
        step = 14 if ndays > 35 else 7 if ndays > 16 else 4
        ticks = pd.date_range(t0, t1, freq=f"{step}D")
        col_ax.xaxis.set_major_locator(
            matplotlib.ticker.FixedLocator(mdates.date2num(ticks))
        )
        # Shared y-range for the column, set by its windiest shared station
        for ax in shared:
            ax.set_ylim(0, max(shared.values()))

    # Fill the empty cells of the column with the most free space with the
    # color wheel and the shared legend; fall back to space below the
    # figure if every cell is used
    short = min(range(ncol), key=lambda i: (len(cols[i]), i))
    n_empty = nrow - len(cols[short])
    if n_empty:
        area = gs[nrow - n_empty:, short].subgridspec(
            1, 2, width_ratios=[1, 1.3], wspace=0.25
        )
        ax_wheel = fig.add_subplot(area[0], projection="polar")
        ax_leg = fig.add_subplot(area[1])
    else:
        ax_wheel = fig.add_axes([0.13, -0.16, 0.12, 0.14],
                                projection="polar")
        ax_leg = fig.add_axes([0.3, -0.16, 0.5, 0.14])
    plot_direction_wheel(ax_wheel, fontsize=fs)
    ax_wheel.text(0.5, -0.2, f"gray: calm\n(<{calm:g} m/s)", fontsize=fs - 1,
                  ha="center", va="top", transform=ax_wheel.transAxes)

    ax_leg.axis("off")
    handles, labels = first_ax.get_legend_handles_labels()
    leg = ax_leg.legend(handles, labels, loc="center left", fontsize=fs,
                        frameon=False, borderaxespad=0, handlelength=1.5)
    # Direction-colored legend marker would show one arbitrary hue; use a
    # neutral, larger marker since the wheel explains the colors
    for handle in leg.legend_handles:
        if hasattr(handle, "set_sizes"):
            handle.set_sizes([20])
            handle.set_array(None)
            handle.set_facecolor("0.4")
            handle.set_alpha(0.8)

    fig.supylabel("Average wind speed (m/s)", fontsize=fs + 2, x=0.085)
    fig.savefig(fid_out, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_all_roses(data, fid_out, calm=0.5, width=7.0, fontsize=8):
    """
    Plot every station's wind rose on one figure with a shared scale.

    The figure is drawn at its printed size (``width`` inches, fitting a
    US letter page), so font sizes are the sizes that appear on the page.

    :type data: dict
    :param data: ``{station: DataFrame}`` of averaged wind data
    :type fid_out: str
    :param fid_out: output figure path
    :type calm: float
    :param calm: calm threshold in m/s
    :type width: float
    :param width: figure width in inches
    :type fontsize: float
    :param fontsize: base font size in points; station names are 2 pt
        larger, radial labels 1.5 pt smaller
    """
    # Order panels by deployment start time, then station code for ties
    stas = sorted(data, key=lambda s: (data[s].index[0], s))
    peaks = [rose_histogram(data[s], calm)[0].sum(axis=0).max()
             for s in stas]
    rmax = np.ceil(max(peaks) * 1.05)
    ncol = min(4, len(stas))
    nrow = int(np.ceil(len(stas) / ncol))
    # Each row is a square rose plus room for the title and calm note
    panel = width / ncol
    height = nrow * panel * 1.2 + 0.4
    fig, axs = plt.subplots(nrow, ncol, figsize=(width, height),
                            subplot_kw={"projection": "polar"},
                            squeeze=False)
    fig.subplots_adjust(left=0.03, right=0.97, top=0.9, bottom=0.11,
                        hspace=0.42, wspace=0.32)
    for ax, sta in zip(axs.flat, stas):
        plot_rose(ax, data[sta], calm=calm, rmax=rmax, legend=False,
                  fontsize=fontsize, calm_y=None)
        # Calm percentage in the lower-left corner, clear of the S label;
        # the threshold itself is given once in the legend title
        pct_calm = rose_histogram(data[sta], calm)[1]
        ax.text(-0.04, -0.04, f"calm\n{pct_calm:.1f}%",
                transform=ax.transAxes, ha="left", va="bottom",
                fontsize=fontsize, linespacing=1.1)
        # Cardinal labels only, so neighbouring E and W labels do not touch;
        # the 45 degree grid lines still mark the intercardinal directions
        ax.set_xticklabels(["N", "", "E", "", "S", "", "W", ""],
                           fontsize=fontsize)
        ax.tick_params(axis="x", pad=-2)
        ax.tick_params(axis="y", labelsize=fontsize - 1.5)
        ax.set_rlabel_position(157.5)
        ax.text(0.5, 1.2, sta, transform=ax.transAxes, ha="center",
                va="bottom", fontsize=fontsize + 2, fontweight="bold")
        ax.text(0.5, 1.13, date_span_short(data[sta]),
                transform=ax.transAxes, ha="center", va="bottom",
                fontsize=fontsize)
    for ax in axs.flat[len(stas):]:
        ax.set_visible(False)
    handles, labels = axs.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels,
               title=f"Speed (m/s); calm is < {calm:g} m/s",
               loc="lower center",
               ncol=len(labels), frameon=False, fontsize=fontsize,
               title_fontsize=fontsize, bbox_to_anchor=(0.5, -0.01),
               handlelength=1.6, columnspacing=1.2)
    fig.savefig(fid_out, dpi=300, bbox_inches="tight", pad_inches=0.03)
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

        daily_stats(df, args.calm).round(3).to_csv(
            os.path.join(args.out, f"{sta}_daily_wind.csv"),
            index_label="date"
        )

    if len(data) > 1:
        plot_all_roses(data, os.path.join(args.out, "ALL_windrose.png"),
                       calm=args.calm)
    if data:
        plot_speed_overview(data, os.path.join(args.out,
                                               "ALL_wind_speed.png"),
                            calm=args.calm, avg_sec=args.avg)
    print(f"figures written to {args.out}")


if __name__ == "__main__":
    main()
