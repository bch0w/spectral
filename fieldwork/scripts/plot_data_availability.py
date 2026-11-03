"""Scan a directory of day-file miniseed data and plot data availability.

Walks every miniseed file in a data directory named in the standard SEED
day-file convention (NETWORK.STATION.LOCATION.CHANNEL.YEAR.DAYOFYEAR),
bins each channel's samples into fixed-width time windows (1 hour, 1
day, etc.), and renders a station/component-by-time 2D grid figure:
cell color encodes what fraction of the window's expected samples are
actually present, and each non-empty cell is annotated with the mean
data value observed in that window. Intended as a "what's available to
use" overview figure for a data-mining manuscript.

Each file is read once in full (fast for day-length miniseed files) and
its samples are binned directly with numpy slicing -- no separate
headonly pass. Per-channel/per-bin statistics can optionally be cached
to a CSV so the figure can be restyled without rescanning the raw data.

Examples:
    python data_avail.py \
        --data-dir /Users/prof/Data/nodes/2025-03_HAVO_VENUS/nodedata \
        --bin-hours 24 --out data_availability.png

    python data_avail.py --data-dir ./nodedata --bin-hours 1 \
        --max-annotated-cells 2000 --out hourly_availability.png \
        --stats-cache cache/nodedata_hourly.csv

    # Re-plot from a cached scan without touching the raw data again
    python scripts/data_avail.py \
        --stats-cache cache/nodedata_daily.csv \
        --out manuscript/figures/havo_data_avail.png \
        --n-color-bins 10 \
        --no-annotate
"""
import argparse
import re
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from obspy import read

# NETWORK.STATION.LOCATION.CHANNEL.YEAR.DAYOFYEAR,
# e.g. "VH.DESD.CV.EHZ.2025.086"
FNAME_RE = re.compile(
    r"^(?P<net>[A-Za-z0-9]+)\.(?P<sta>[A-Za-z0-9]+)\."
    r"(?P<loc>[A-Za-z0-9]*)\.(?P<chan>[A-Za-z0-9]+)\."
    r"(?P<year>\d{4})\.(?P<doy>\d{3})$"
)


def discover_files(data_dir, pattern):
    """
    Find and parse every miniseed day-file in `data_dir` matching
    `FNAME_RE`. Files that don't match the naming convention are
    skipped with a warning rather than raising, since a data directory
    can legitimately contain sidecar files (logs, checksums, etc).

    :type data_dir: str or pathlib.Path
    :param data_dir: directory to scan (non-recursive)
    :type pattern: str
    :param pattern: glob pattern used to list candidate files, e.g. "*"

    :rtype: pandas.DataFrame
    :return: one row per file, columns seed_id, net, sta, loc, chan,
        year, doy, path
    """
    manifest = []
    for path in sorted(Path(data_dir).glob(pattern)):
        if not path.is_file():
            continue
        match = FNAME_RE.match(path.name)
        if not match:
            print(f"skipping {path.name}: not NET.STA.LOC.CHAN.YYYY.DDD")
            continue
        seed_id = (f"{match['net']}.{match['sta']}."
                   f"{match['loc']}.{match['chan']}")
        manifest.append({
            "seed_id": seed_id, "net": match["net"], "sta": match["sta"],
            "loc": match["loc"], "chan": match["chan"],
            "year": int(match["year"]), "doy": int(match["doy"]),
            "path": path,
        })
    if not manifest:
        raise FileNotFoundError(
            "No miniseed files matching NET.STA.LOC.CHAN.YYYY.DDD found "
            f"in {data_dir}")
    return pd.DataFrame(manifest)


def build_bin_edges(manifest, bin_hours):
    """
    Build a global, UTC-day-aligned array of bin edges spanning every
    file in `manifest`, from the first file's day (00:00:00) through
    the day after the last file's day.

    :type manifest: pandas.DataFrame
    :param manifest: output of `discover_files`
    :type bin_hours: float
    :param bin_hours: bin width in hours

    :rtype: numpy.ndarray
    :return: 1D array of Unix timestamps (float seconds), length
        n_bins + 1
    """
    from obspy import UTCDateTime

    year_doy = manifest[["year", "doy"]].apply(tuple, axis=1)
    start_row = manifest.loc[year_doy.idxmin()]
    end_row = manifest.loc[year_doy.idxmax()]
    global_start = UTCDateTime(year=int(start_row["year"]),
                                julday=int(start_row["doy"]))
    global_end = UTCDateTime(year=int(end_row["year"]),
                              julday=int(end_row["doy"])) + 86400

    bin_seconds = bin_hours * 3600
    span = global_end.timestamp - global_start.timestamp
    n_bins = int(np.ceil(span / bin_seconds))
    return global_start.timestamp + np.arange(n_bins + 1) * bin_seconds


def scan_availability(manifest, bin_edges):
    """
    Read every file in `manifest` once and accumulate per-channel,
    per-bin sample counts and min/sum/max statistics.

    :type manifest: pandas.DataFrame
    :param manifest: output of `discover_files`
    :type bin_edges: numpy.ndarray
    :param bin_edges: output of `build_bin_edges`

    :rtype: tuple
    :return: (seed_ids, samples_present, expected_samples, min_val,
        mean_val, max_val) -- the first is a sorted list of channel
        ids, the rest are (n_channels, n_bins) numpy arrays
    """
    seed_ids = sorted(manifest["seed_id"].unique())
    row_of = {seed_id: i for i, seed_id in enumerate(seed_ids)}
    n_rows, n_bins = len(seed_ids), len(bin_edges) - 1

    samples_present = np.zeros((n_rows, n_bins), dtype=np.int64)
    sum_val = np.zeros((n_rows, n_bins), dtype=np.float64)
    min_val = np.full((n_rows, n_bins), np.inf)
    max_val = np.full((n_rows, n_bins), -np.inf)
    sampling_rate = np.full(n_rows, np.nan)

    n_files = len(manifest)
    for k, file_row in enumerate(manifest.itertuples(), start=1):
        print(f"[{k}/{n_files}] {file_row.path.name}")
        stream = read(str(file_row.path))
        r = row_of[file_row.seed_id]

        for trace in stream:
            sr = trace.stats.sampling_rate
            if np.isnan(sampling_rate[r]):
                sampling_rate[r] = sr
            elif not np.isclose(sampling_rate[r], sr):
                print(f"  warning: {file_row.seed_id} sampling rate "
                      f"changed {sampling_rate[r]} -> {sr} Hz in "
                      f"{file_row.path.name}")

            t0 = trace.stats.starttime.timestamp
            t1 = trace.stats.endtime.timestamp + trace.stats.delta
            data = trace.data.astype(np.float64, copy=False)

            i_start = np.searchsorted(bin_edges, t0, side="right") - 1
            i_start = max(0, int(i_start))
            i_end = min(n_bins, int(np.searchsorted(bin_edges, t1)))
            for i in range(i_start, i_end):
                seg_start = max(bin_edges[i], t0)
                seg_end = min(bin_edges[i + 1], t1)
                if seg_end <= seg_start:
                    continue
                idx0 = max(0, int(round((seg_start - t0) * sr)))
                idx1 = min(data.size, int(round((seg_end - t0) * sr)))
                if idx1 <= idx0:
                    continue
                seg = data[idx0:idx1]
                samples_present[r, i] += seg.size
                sum_val[r, i] += seg.sum()
                min_val[r, i] = min(min_val[r, i], seg.min())
                max_val[r, i] = max(max_val[r, i], seg.max())
        del stream

    bin_seconds = bin_edges[1] - bin_edges[0]
    expected_samples = np.outer(sampling_rate,
                                 np.full(n_bins, bin_seconds))
    with np.errstate(invalid="ignore"):
        safe_n = np.maximum(samples_present, 1)
        mean_val = np.where(samples_present > 0, sum_val / safe_n, np.nan)
    min_val = np.where(samples_present > 0, min_val, np.nan)
    max_val = np.where(samples_present > 0, max_val, np.nan)

    return (seed_ids, samples_present, expected_samples, min_val,
            mean_val, max_val)


def to_tidy_dataframe(seed_ids, bin_edges, fraction, min_val, mean_val,
                       max_val):
    """
    Flatten the (n_channels, n_bins) arrays into a long DataFrame for
    CSV caching.
    """
    n_rows, n_bins = fraction.shape
    bin_start = pd.to_datetime(bin_edges[:-1], unit="s")
    bin_end = pd.to_datetime(bin_edges[1:], unit="s")
    return pd.DataFrame({
        "seed_id": np.repeat(seed_ids, n_bins),
        "bin_start": np.tile(bin_start, n_rows),
        "bin_end": np.tile(bin_end, n_rows),
        "fraction": fraction.ravel(),
        "min": min_val.ravel(),
        "mean": mean_val.ravel(),
        "max": max_val.ravel(),
    })


def from_tidy_dataframe(df):
    """
    Reconstruct the (n_channels, n_bins) arrays a cached CSV produced
    by `to_tidy_dataframe`.
    """
    seed_ids = sorted(df["seed_id"].unique())
    bin_start = np.sort(df["bin_start"].unique())
    bin_end = np.sort(df["bin_end"].unique())
    edges = np.append(bin_start, bin_end[-1])
    bin_edges = edges.astype("datetime64[s]").astype("int64")
    bin_edges = bin_edges.astype(np.float64)

    def _pivot(col):
        table = df.pivot(index="seed_id", columns="bin_start", values=col)
        return table.loc[seed_ids].to_numpy()

    fraction = _pivot("fraction")
    min_val = _pivot("min")
    mean_val = _pivot("mean")
    max_val = _pivot("max")
    return seed_ids, bin_edges, fraction, min_val, mean_val, max_val


def condense_channels(seed_ids, fraction):
    """
    Collapse the component channels of each instrument into a single
    row by replacing the last character of the channel code with '?',
    e.g. EHZ/EHN/EHE -> EH? and LWD/LWS -> LW?. The condensed row takes
    the minimum availability across its components in each bin, so a
    gap on any one component is still visible.

    :type seed_ids: list of str
    :param seed_ids: sorted NET.STA.LOC.CHAN row labels
    :type fraction: numpy.ndarray
    :param fraction: (n_channels, n_bins) fraction of expected samples
        present

    :rtype: tuple
    :return: (condensed_ids, condensed_fraction) -- a sorted list of
        NET.STA.LOC.CH? labels and the matching (n_instruments, n_bins)
        array
    """
    groups = {}
    for i, seed_id in enumerate(seed_ids):
        groups.setdefault(f"{seed_id[:-1]}?", []).append(i)
    condensed_ids = sorted(groups)
    condensed = np.array([np.nanmin(fraction[groups[gid]], axis=0)
                          for gid in condensed_ids])
    return condensed_ids, condensed


def trim_time_range(bin_edges, arrays, starttime=None, endtime=None):
    """
    Restrict the bins to a time window, keeping every bin that overlaps
    [starttime, endtime). Trimming the data (rather than only setting
    the axis limits) keeps the colorscale, ticks and title consistent
    with what is shown.

    :type bin_edges: numpy.ndarray
    :param bin_edges: 1D array of Unix timestamps, length n_bins + 1
    :type arrays: list of numpy.ndarray
    :param arrays: (n_channels, n_bins) arrays to trim along the bin axis
    :type starttime: str or None
    :param starttime: UTC start of the window, any format understood by
        pandas.Timestamp, e.g. "2025-03-25" or "2025-03-25T12:00"
    :type endtime: str or None
    :param endtime: UTC end of the window (exclusive), same format

    :rtype: tuple
    :return: (trimmed_bin_edges, trimmed_arrays)
    """
    t0 = pd.Timestamp(starttime).timestamp() if starttime else -np.inf
    t1 = pd.Timestamp(endtime).timestamp() if endtime else np.inf
    keep = (bin_edges[1:] > t0) & (bin_edges[:-1] < t1)
    if not keep.any():
        raise ValueError(
            f"No bins between --starttime {starttime} and --endtime "
            f"{endtime}; data spans "
            f"{pd.to_datetime(bin_edges[0], unit='s')} to "
            f"{pd.to_datetime(bin_edges[-1], unit='s')}")
    idx = np.flatnonzero(keep)
    trimmed_edges = bin_edges[idx[0]:idx[-1] + 2]
    return trimmed_edges, [arr[:, keep] for arr in arrays]


def month_restart_ticks(bin_edges, interval_days):
    """
    Build tick positions every `interval_days` days, restarting the
    count at each calendar month boundary -- the month's 1st is day
    index 0, so the next tick after it falls on day `interval_days`
    (index `interval_days - 1`), not a continuation of the previous
    month's spacing. The 1st of every month is always included as a
    tick regardless of `interval_days`, so month boundaries stay
    visible. E.g. for interval_days=7 in April: Apr 01, 07, 14, 21,
    28, then May resets and starts again at May 01, 07.

    :type bin_edges: numpy.ndarray
    :param bin_edges: 1D array of Unix timestamps
    :type interval_days: int
    :param interval_days: spacing between ticks within a month, in days

    :rtype: pandas.DatetimeIndex
    :return: sorted tick positions
    """
    start = pd.to_datetime(bin_edges[0], unit="s").normalize()
    end = pd.to_datetime(bin_edges[-1], unit="s")
    month_start = start.replace(day=1)

    ticks = []
    while month_start <= end:
        next_month = month_start + pd.offsets.MonthBegin(1)
        first_tick = month_start + pd.Timedelta(days=interval_days - 1)
        ticks.append(month_start)
        ticks.extend(pd.date_range(first_tick, next_month,
                                    freq=f"{interval_days}D",
                                    inclusive="left"))
        month_start = next_month
    return pd.DatetimeIndex(sorted(set(ticks)))


def julian_day_range_label(bin_edges):
    """
    Describe the Julian-day (day-of-year) span covered by `bin_edges`,
    for use in the figure title.

    :type bin_edges: numpy.ndarray
    :param bin_edges: 1D array of Unix timestamps, length n_bins + 1

    :rtype: str
    :return: e.g. "DOY 082-123, 2025", or "DOY 360, 2025 - DOY 010,
        2026" if the span crosses a year boundary
    """
    start_dt = pd.to_datetime(bin_edges[0], unit="s")
    end_dt = pd.to_datetime(bin_edges[-1], unit="s") - pd.Timedelta(seconds=1)
    if start_dt.year == end_dt.year:
        return (f"DOY {start_dt.dayofyear:03d}-{end_dt.dayofyear:03d}, "
                f"{start_dt.year}")
    return (f"DOY {start_dt.dayofyear:03d}, {start_dt.year} - "
            f"DOY {end_dt.dayofyear:03d}, {end_dt.year}")


def plot_availability(seed_ids, bin_edges, fraction, min_val, mean_val,
                       max_val, bin_hours, out, annotate=True, title=None,
                       max_annotated_cells=800, annotate_fontsize=6.5,
                       n_color_bins=11, tick_interval_days=5, dpi=200,
                       segments=None):
    """
    Render the station/component-by-time availability grid and save
    it. If `segments` is given, each time segment gets its own panel,
    with broken-axis marks between them.

    :type seed_ids: list of str
    :param seed_ids: sorted NET.STA.LOC.CHAN row labels
    :type bin_edges: numpy.ndarray
    :param bin_edges: 1D array of Unix timestamps, length n_bins + 1
    :type fraction: numpy.ndarray
    :param fraction: (n_channels, n_bins) fraction of expected samples
        present
    :type min_val: numpy.ndarray
    :type mean_val: numpy.ndarray
    :type max_val: numpy.ndarray
    :param min_val, mean_val, max_val: (n_channels, n_bins) per-bin
        data value statistics, NaN where `fraction` is 0
    :type bin_hours: float
    :param bin_hours: bin width in hours, used only for the title
    :type out: str or pathlib.Path
    :param out: image path to save the figure to
    :type annotate: bool
    :param annotate: annotate non-empty cells with their mean value
    :type max_annotated_cells: int
    :param max_annotated_cells: skip annotation (with a warning) if
        the grid has more cells than this -- text stops being legible
        well before it stops being renderable
    :type annotate_fontsize: float
    :param annotate_fontsize: font size for the mean-value cell text
    :type n_color_bins: int
    :param n_color_bins: number of discrete color levels the
        availability colorscale is split into, spanning the minimum
        non-zero fraction present up to 1 (default: 11). Complete
        bins (fraction == 1) are drawn in a separate color instead
    :type tick_interval_days: int
    :param tick_interval_days: spacing, in days, between x-axis tick
        marks within a calendar month (see `month_restart_ticks`)
    :type dpi: int
    :param dpi: figure save resolution
    :type segments: list of tuple or None
    :param segments: [(first_bin, last_bin + 1), ...] index ranges to
        plot as separate panels (see `find_segments`); None plots every
        bin in one panel
    """
    n_rows, n_bins = fraction.shape
    if segments is None:
        segments = [(0, n_bins)]
    seg_bins = [i1 - i0 for i0, i1 in segments]

    fig_w = float(np.clip(2 + sum(seg_bins) * 0.3, 8, 40))
    fig_h = float(np.clip(1.5 + n_rows * 0.35, 4, 20))
    # One panel per segment, sized by its number of bins, sharing rows
    fig, axes = plt.subplots(
        1, len(segments), figsize=(fig_w, fig_h), sharey=True,
        squeeze=False, layout="constrained",
        gridspec_kw={"width_ratios": seg_bins})
    fig.get_layout_engine().set(wspace=0.02, w_pad=0.02)
    axes = axes[0]

    cmap = plt.get_cmap("RdYlGn", n_color_bins).copy()
    cmap.set_bad(color="#d9d9d9")
    # Start the colorscale at the lowest non-zero availability in the set
    # (floored to 2 decimals for tidy ticks) rather than at 0, so the
    # color levels resolve the range of availability actually present.
    # Empty bins are masked and drawn gray, so they don't need a color.
    present = fraction[fraction > 0]
    vmin = np.floor(present.min() * 100) / 100 if present.size else 0.
    if vmin >= 1:
        vmin = 0.
    bounds = np.linspace(vmin, 1, n_color_bins + 1)
    # Use just enough decimals that adjacent colorbar ticks stay distinct
    # when the range is narrow (e.g. 0.99-1.00 needs 4 decimals)
    step = bounds[1] - bounds[0]
    n_decimals = max(2, int(np.ceil(-np.log10(step))))
    norm = mcolors.BoundaryNorm(bounds, n_color_bins)
    complete_color = "#2166ac"

    total_cells = n_rows * sum(seg_bins)
    do_annotate = annotate and total_cells <= max_annotated_cells
    if annotate and not do_annotate:
        print(f"Skipping cell annotations: grid has {total_cells} "
              f"cells > --max-annotated-cells {max_annotated_cells} "
              "(coarsen --bin-hours, or raise --max-annotated-cells, "
              "to force annotation).")

    for k, (ax, (i0, i1)) in enumerate(zip(axes, segments)):
        # Only the last panel labels its end date; elsewhere it would
        # collide with the next panel's start date across the break
        mesh = draw_segment(ax, seed_ids, bin_edges[i0:i1 + 1],
                            fraction[:, i0:i1], mean_val[:, i0:i1], cmap,
                            norm, complete_color, do_annotate,
                            annotate_fontsize, tick_interval_days,
                            end_tick=k == len(segments) - 1)

    axes[0].set_yticks(np.arange(n_rows) + 0.5)
    axes[0].set_yticklabels(seed_ids)
    axes[0].invert_yaxis()  # shared, so this flips every panel
    if len(axes) > 1:
        draw_axis_breaks(axes)
    axes[-1].legend(handles=[mpatches.Patch(color=complete_color,
                                            label="Complete (100%)")],
                    loc="lower right", bbox_to_anchor=(1.0, 1.0),
                    fontsize=10, frameon=False, borderaxespad=0.5)

    fig.supxlabel("Time (UTC)")
    day_range = julian_day_range_label(bin_edges)
    title = title or (f"Data Availability ({bin_hours:g}-hour bins, "
                      f"{day_range})")
    fig.suptitle(title)

    cbar = fig.colorbar(mesh, ax=list(axes), pad=0.01, boundaries=bounds,
                         ticks=bounds)
    # Top class is [x, 1) -- label its upper edge "< 1" so it isn't
    # mistaken for the separately colored complete (100%) bins
    tick_labels = [f"{b:.{n_decimals}f}" for b in bounds[:-1]]
    cbar.set_ticklabels(tick_labels + ["< 1"])
    cbar.set_label("Fraction of Expected Samples")

    if do_annotate:
        fig.text(0.01, 0.01,
                "Cell text: mean data value in that bin (raw miniseed "
                "units)",
                fontsize=8, color="dimgray")
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"Saved figure to {out}")


def draw_segment(ax, seed_ids, bin_edges, fraction, mean_val, cmap, norm,
                 complete_color, annotate, annotate_fontsize,
                 tick_interval_days, end_tick=True):
    """
    Draw one continuous time segment of the availability grid on `ax`.

    :type ax: matplotlib.axes.Axes
    :param ax: axes to draw on
    :type seed_ids: list of str
    :param seed_ids: row labels, used for the station separator lines
    :type bin_edges: numpy.ndarray
    :param bin_edges: this segment's Unix timestamp bin edges
    :type fraction: numpy.ndarray
    :param fraction: (n_channels, n_seg_bins) fraction of expected
        samples present
    :type mean_val: numpy.ndarray
    :param mean_val: (n_channels, n_seg_bins) mean data value per bin
    :type cmap: matplotlib.colors.Colormap
    :param cmap: colormap for partially complete bins
    :type norm: matplotlib.colors.BoundaryNorm
    :param norm: color normalization shared by all segments
    :type complete_color: str
    :param complete_color: color for complete (100%) bins
    :type annotate: bool
    :param annotate: annotate non-empty cells with their mean value
    :type annotate_fontsize: float
    :param annotate_fontsize: font size for the mean-value cell text
    :type tick_interval_days: int
    :param tick_interval_days: x tick spacing in days within a month
    :type end_tick: bool
    :param end_tick: label the segment's end edge with a tick
    :rtype: matplotlib.collections.QuadMesh
    :return: the partial-availability mesh, for the colorbar
    """
    n_rows, n_bins = fraction.shape
    x_edges = pd.to_datetime(bin_edges, unit="s")
    y_edges = np.arange(n_rows + 1)

    # Bins with every expected sample present get their own color, drawn
    # as a separate layer, so even a few missing samples (which would
    # otherwise land in the top color class) are visibly incomplete
    complete = fraction >= 1
    frac_masked = np.ma.masked_where((fraction <= 0) | complete, fraction)
    mesh = ax.pcolormesh(x_edges, y_edges, frac_masked, cmap=cmap,
                          norm=norm, edgecolors="white", linewidth=0.4)
    ax.pcolormesh(x_edges, y_edges,
                  np.ma.masked_where(~complete, np.ones_like(fraction)),
                  cmap=mcolors.ListedColormap([complete_color]),
                  edgecolors="white", linewidth=0.4)

    if annotate:
        centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_centers = pd.to_datetime(centers, unit="s")
        for r in range(n_rows):
            for c in range(n_bins):
                if fraction[r, c] <= 0:
                    continue
                text = f"{mean_val[r, c]:.3g}"
                ax.text(bin_centers[c], r + 0.5, text, ha="center",
                        va="center", fontsize=annotate_fontsize,
                        color="black",
                        bbox=dict(facecolor="white", alpha=0.55,
                                  edgecolor="none", pad=0.5))

    # Separator lines between rows: heavier where the station changes,
    # lighter where just the location code changes (e.g. DESD.EX ->
    # DESD.WS) -- never where only the component/channel changes.
    stations = [seed_id.split(".")[1] for seed_id in seed_ids]
    locations = [seed_id.split(".")[2] for seed_id in seed_ids]
    for r in range(1, n_rows):
        if stations[r] != stations[r - 1]:
            ax.axhline(r, color="black", linewidth=1.5)
        elif locations[r] != locations[r - 1]:
            ax.axhline(r, color="black", linewidth=0.8)

    ticks = month_restart_ticks(bin_edges, tick_interval_days)
    last = x_edges[-1] if end_tick else x_edges[-1] - pd.Timedelta(days=1)
    ticks = ticks[(ticks >= x_edges[0]) & (ticks <= last)]
    range_ticks = [x_edges[0], x_edges[-1]] if end_tick else [x_edges[0]]
    ticks = pd.DatetimeIndex(sorted(set(ticks) | set(range_ticks)))
    ax.set_xticks(ticks)
    tick_labels = [f"{t.strftime('%b %d')}\n{t.dayofyear:03d}"
                   for t in ticks]
    ax.set_xticklabels(tick_labels, rotation=45, ha="center")
    ax.set_xlim(x_edges[0], x_edges[-1])
    return mesh


def draw_axis_breaks(axes, size=10):
    """
    Join side-by-side panels with broken-axis marks: hide the facing
    spines and draw a pair of diagonal slashes at the top and bottom of
    each break.

    :type axes: list of matplotlib.axes.Axes
    :param axes: panels ordered left to right
    :type size: float
    :param size: slash marker size in points
    """
    slash = dict(marker=[(-1, -1), (1, 1)], markersize=size,
                 linestyle="none", color="black", mec="black", mew=1.2,
                 clip_on=False)
    for left, right in zip(axes[:-1], axes[1:]):
        left.spines["right"].set_visible(False)
        right.spines["left"].set_visible(False)
        right.tick_params(axis="y", left=False)
        left.plot([1, 1], [0, 1], transform=left.transAxes, **slash)
        right.plot([0, 0], [0, 1], transform=right.transAxes, **slash)


def sort_rows(seed_ids, arrays, sort_by="name", reverse=False):
    """
    Reorder the rows by name, or by when each row's data starts or ends.
    Ties (e.g. the components of one instrument, which usually share
    times) fall back to name order so related rows stay together, and
    rows with no data at all always go last.

    :type seed_ids: list of str
    :param seed_ids: row labels, in name order
    :type arrays: list of numpy.ndarray
    :param arrays: (n_rows, n_bins) arrays to reorder along the row axis;
        the first must be the availability fraction
    :type sort_by: str
    :param sort_by: "name", "starttime" (first bin with data) or
        "endtime" (last bin with data)
    :type reverse: bool
    :param reverse: sort latest (or Z-A) first instead
    :rtype: tuple
    :return: (sorted_seed_ids, sorted_arrays)
    """
    has_data = np.nan_to_num(arrays[0]) > 0
    n_bins = has_data.shape[1]
    keys = []
    for r, seed_id in enumerate(seed_ids):
        bins = np.flatnonzero(has_data[r])
        if sort_by == "starttime":
            t = bins[0] if bins.size else None
        elif sort_by == "endtime":
            t = bins[-1] if bins.size else None
        else:
            t = 0
        if t is not None and reverse and sort_by != "name":
            t = n_bins - t  # flip the time order, keep names A-Z
        keys.append((t is None, t or 0, seed_id))
    order = sorted(range(len(seed_ids)), key=lambda r: keys[r])
    if reverse and sort_by == "name":
        order = order[::-1]
    return ([seed_ids[r] for r in order],
            [arr[order] for arr in arrays])


def find_segments(fraction, bin_edges, min_gap_days):
    """
    Split the bins into continuous segments, cutting out every run of
    bins where no row has any data that lasts at least `min_gap_days`.
    Shorter empty runs stay in place (drawn gray) so brief outages are
    still visible.

    :type fraction: numpy.ndarray
    :param fraction: (n_channels, n_bins) fraction of expected samples
        present
    :type bin_edges: numpy.ndarray
    :param bin_edges: 1D array of Unix timestamps, length n_bins + 1
    :type min_gap_days: float
    :param min_gap_days: shortest all-empty span, in days, to cut out
    :rtype: list of tuple
    :return: [(first_bin, last_bin + 1), ...] index ranges to plot
    """
    empty = ~(np.nan_to_num(fraction) > 0).any(axis=0)
    segments = []
    start = 0
    i = 0
    n_bins = empty.size
    while i < n_bins:
        if not empty[i]:
            i += 1
            continue
        j = i
        while j < n_bins and empty[j]:
            j += 1
        gap_days = (bin_edges[j] - bin_edges[i]) / 86400
        if gap_days >= min_gap_days:
            if i > start:
                segments.append((start, i))
            print(f"Cutting {gap_days:g}-day gap with no data: "
                  f"{pd.to_datetime(bin_edges[i], unit='s')} to "
                  f"{pd.to_datetime(bin_edges[j], unit='s')}")
            start = j
        i = j
    if start < n_bins:
        segments.append((start, n_bins))
    if not segments:
        raise ValueError("No data in any bin to plot")
    return segments


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--data-dir", default=None,
        help="directory of NET.STA.LOC.CHAN.YYYY.DDD miniseed day-files")
    parser.add_argument(
        "--pattern", default="*",
        help="glob pattern for candidate files in --data-dir")
    parser.add_argument(
        "--bin-hours", type=float, default=24.0,
        help="width of each time bin in hours, e.g. 1 or 24 (default: 24)")
    parser.add_argument("--out", default="data_availability.png",
                         help="output figure path")
    parser.add_argument("--dpi", type=int, default=200)
    parser.add_argument(
        "--annotate", action=argparse.BooleanOptionalAction, default=True,
        help="annotate non-empty cells with their mean value "
             "(default: on)")
    parser.add_argument(
        "--max-annotated-cells", type=int, default=800,
        help="skip annotation if the grid exceeds this many cells "
             "(default: 800)")
    parser.add_argument("--annotate-fontsize", type=float, default=6.5)
    parser.add_argument(
        "--n-color-bins", type=int, default=11,
        help="number of discrete color levels for the availability "
             "colorscale, from the minimum non-zero fraction to 1 "
             "(default: 11)")
    parser.add_argument("--title", type=str, default=None)
    parser.add_argument(
        "--starttime", default=None,
        help="UTC start of the plotted window, e.g. 2025-03-25 or "
             "2025-03-25T12:00; earlier data are not shown")
    parser.add_argument(
        "--endtime", default=None,
        help="UTC end of the plotted window (exclusive), e.g. "
             "2025-04-10 shows data through the end of Apr 09")
    parser.add_argument(
        "--condense", action="store_true",
        help="show one row per instrument (e.g. EH? instead of "
             "EHZ/EHN/EHE), taking the minimum availability across its "
             "components; disables cell annotation since mean values "
             "aren't comparable across components")
    parser.add_argument(
        "--sort", choices=["name", "starttime", "endtime"], default="name",
        help="row order: alphabetical by name (default), or by the first "
             "(starttime) or last (endtime) bin with data; ties are "
             "ordered by name and rows with no data go last")
    parser.add_argument(
        "--reverse-sort", action="store_true",
        help="reverse --sort, e.g. latest starttime first")
    parser.add_argument(
        "--break-gaps", type=float, nargs="?", const=7.0, default=None,
        metavar="DAYS",
        help="cut out spans with no data on any row lasting at least "
             "DAYS (default 7 if no value is given) and draw the "
             "remaining periods as separate panels with a broken axis")
    parser.add_argument(
        "--tick-interval-days", type=int, default=7,
        help="x-axis tick spacing in days within a calendar month, "
             "restarting at each month's 1st (default: 5)")
    parser.add_argument(
        "--stats-cache", default=None,
        help="CSV path to save (or load, if it already exists) "
             "per-bin statistics, so the figure can be restyled "
             "without rescanning")
    parser.add_argument(
        "--refresh", action="store_true",
        help="rescan the raw data even if --stats-cache already exists")
    args = parser.parse_args()

    have_cache = (args.stats_cache and Path(args.stats_cache).exists()
                  and not args.refresh)
    if args.data_dir is None and not have_cache:
        parser.error(
            "--data-dir is required unless loading an existing "
            "--stats-cache")
    return args


def main():
    args = parse_args()

    cache_path = Path(args.stats_cache) if args.stats_cache else None
    if cache_path and cache_path.exists() and not args.refresh:
        print(f"Loading cached statistics from {cache_path}")
        df = pd.read_csv(cache_path, parse_dates=["bin_start", "bin_end"])
        seed_ids, bin_edges, fraction, min_val, mean_val, max_val = \
            from_tidy_dataframe(df)
        span = df["bin_end"].iloc[0] - df["bin_start"].iloc[0]
        bin_hours = span.total_seconds() / 3600
    else:
        manifest = discover_files(args.data_dir, args.pattern)
        print(f"Found {len(manifest)} files across "
              f"{manifest['seed_id'].nunique()} channels")

        bin_edges = build_bin_edges(manifest, args.bin_hours)
        (seed_ids, samples_present, expected_samples, min_val, mean_val,
         max_val) = scan_availability(manifest, bin_edges)
        with np.errstate(invalid="ignore"):
            fraction = np.clip(samples_present / expected_samples, 0, 1)
        bin_hours = args.bin_hours

        if cache_path:
            df = to_tidy_dataframe(seed_ids, bin_edges, fraction, min_val,
                                    mean_val, max_val)
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(cache_path, index=False)
            print(f"Saved per-bin statistics to {cache_path}")

    if args.starttime or args.endtime:
        bin_edges, (fraction, min_val, mean_val, max_val) = \
            trim_time_range(bin_edges,
                            [fraction, min_val, mean_val, max_val],
                            args.starttime, args.endtime)

    annotate = args.annotate
    if args.condense:
        n_channels = len(seed_ids)
        seed_ids, fraction = condense_channels(seed_ids, fraction)
        print(f"Condensed {n_channels} channels into {len(seed_ids)} "
              "instruments")
        if annotate:
            print("Skipping cell annotations: mean values aren't "
                  "comparable across condensed components")
            annotate = False
        min_val = mean_val = max_val = np.full_like(fraction, np.nan)

    if args.sort != "name" or args.reverse_sort:
        seed_ids, (fraction, min_val, mean_val, max_val) = sort_rows(
            seed_ids, [fraction, min_val, mean_val, max_val], args.sort,
            args.reverse_sort)

    segments = None
    if args.break_gaps is not None:
        segments = find_segments(fraction, bin_edges, args.break_gaps)

    plot_availability(seed_ids, bin_edges, fraction, min_val, mean_val,
                       max_val, bin_hours, args.out, annotate=annotate,
                       title=args.title,
                       max_annotated_cells=args.max_annotated_cells,
                       annotate_fontsize=args.annotate_fontsize,
                       n_color_bins=args.n_color_bins,
                       tick_interval_days=args.tick_interval_days,
                       dpi=args.dpi, segments=segments)


if __name__ == "__main__":
    main()
