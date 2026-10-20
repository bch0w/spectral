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
                       max_val, bin_hours, out, annotate=True,
                       max_annotated_cells=800, annotate_fontsize=6.5,
                       n_color_bins=11, tick_interval_days=5, dpi=200):
    """
    Render the station/component-by-time availability grid and save
    it.

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
    :param n_color_bins: number of discrete color levels the 0-1
        availability colorscale is split into (default: 11)
    :type tick_interval_days: int
    :param tick_interval_days: spacing, in days, between x-axis tick
        marks within a calendar month (see `month_restart_ticks`)
    :type dpi: int
    :param dpi: figure save resolution
    """
    n_rows, n_bins = fraction.shape
    x_edges = pd.to_datetime(bin_edges, unit="s")
    y_edges = np.arange(n_rows + 1)

    fig_w = float(np.clip(2 + n_bins * 0.3, 8, 40))
    fig_h = float(np.clip(1.5 + n_rows * 0.35, 4, 20))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    cmap = plt.get_cmap("RdYlGn", n_color_bins).copy()
    cmap.set_bad(color="#d9d9d9")
    bounds = np.linspace(0, 1, n_color_bins + 1)
    norm = mcolors.BoundaryNorm(bounds, n_color_bins)
    frac_masked = np.ma.masked_where(fraction <= 0, fraction)
    mesh = ax.pcolormesh(x_edges, y_edges, frac_masked, cmap=cmap,
                          norm=norm, edgecolors="white", linewidth=0.4)

    total_cells = n_rows * n_bins
    do_annotate = annotate and total_cells <= max_annotated_cells
    if annotate and not do_annotate:
        print(f"Skipping cell annotations: grid has {total_cells} "
              f"cells > --max-annotated-cells {max_annotated_cells} "
              "(coarsen --bin-hours, or raise --max-annotated-cells, "
              "to force annotation).")
    if do_annotate:
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

    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels(seed_ids)
    ax.invert_yaxis()
    ticks = month_restart_ticks(bin_edges, tick_interval_days)
    range_ticks = pd.DatetimeIndex([x_edges[0], x_edges[-1]])
    ticks = pd.DatetimeIndex(sorted(set(ticks) | set(range_ticks)))
    ax.set_xticks(ticks)
    tick_labels = [f"{t.strftime('%b %d')}\n{t.dayofyear:03d}"
                   for t in ticks]
    ax.set_xticklabels(tick_labels, rotation=45, ha="right")
    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_xlabel("Time (UTC)")
    day_range = julian_day_range_label(bin_edges)
    ax.set_title(f"Data Availability ({bin_hours:g}-hour "
                 f"bins, {day_range})\n"
                 f"CV=covered, EX=exposed, WS=windshield")

    cbar = fig.colorbar(mesh, ax=ax, pad=0.01, boundaries=bounds,
                         ticks=bounds)
    cbar.set_label("Fraction of expected samples present")

    if do_annotate:
        fig.text(0.01, 0.01,
                "Cell text: mean data value in that bin (raw miniseed "
                "units)",
                fontsize=8, color="dimgray")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"Saved figure to {out}")


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
        help="number of discrete color levels for the 0-1 availability "
             "colorscale (default: 11)")
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

    plot_availability(seed_ids, bin_edges, fraction, min_val, mean_val,
                       max_val, bin_hours, args.out, annotate=args.annotate,
                       max_annotated_cells=args.max_annotated_cells,
                       annotate_fontsize=args.annotate_fontsize,
                       n_color_bins=args.n_color_bins,
                       tick_interval_days=args.tick_interval_days,
                       dpi=args.dpi)


if __name__ == "__main__":
    main()
