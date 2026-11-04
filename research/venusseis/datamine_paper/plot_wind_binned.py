#!/usr/bin/env python
"""
Plot how seismic noise at each station changes with wind speed.

For one site, every station's hourly PSDs are grouped by the mean wind
speed measured during that hour at the site's anemometer (the XX station
with suffix A, e.g. MITDA). For each wind-speed bin the median spectrum is
taken, and each panel shows its difference from the same station's median
spectrum in the calmest bin (0-1 m/s by default), in the style of Seppi et
al. (2026, Figure 9). Rows are stations (Permanent, Covered, Uncovered,
Windscreen) and columns are components (Z, N, E).

PSDs are smoothed onto 1/12-octave bins as in plot_permanent_relative.py
and limited to the site's common time window. Wind is read with
plot_wind.py (faulty samples removed), averaged in 10-minute windows, and
each PSD hour gets the mean of the 10-minute windows it spans.

Example::

    python plot_wind_binned.py --site MITD
    python plot_wind_binned.py --site PAX --bins 0,1,2,4,6,8
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from plot_permanent_relative import (  # noqa: E402
    KEY_ORDER, MUTED_TEXT, channel_data, common_window, key_label,
    load_stations, log_bins, station_key, style_axis, to_bins,
)
from plot_wind import find_wind_files, load_station  # noqa: E402


# Calm reference bin in neutral gray; windier bins on a light -> dark
# warm sequential ramp so color order follows wind speed
CALM_COLOR = "#8a8984"
WIND_CMAP = "YlOrRd"
WIND_AVG_SEC = 600
PSD_WINDOW_SEC = 3600


def hourly_wind(site, data_dir, suffix="A"):
    """
    Load a site's wind speed as 10-minute means keyed by start time.

    :type site: str
    :param site: site code, e.g. MITD
    :type data_dir: str
    :param data_dir: MiniSEED directory holding LWS/LWD day files
    :type suffix: str
    :param suffix: suffix of the anemometer station code
    :rtype: dict
    :return: Unix start time (s) -> mean wind speed (m/s)
    """
    station = f"{site}{suffix}"
    files = find_wind_files(data_dir, [station]).get(station)
    if not files or not files["LWS"] or not files["LWD"]:
        raise SystemExit(f"No wind data for {station} in {data_dir}")
    df = load_station(files["LWS"], files["LWD"], avg_sec=WIND_AVG_SEC)
    times = df.index.asi8 // 1_000_000_000
    return dict(zip(times.tolist(), df.speed.to_numpy()))


def wind_for_segments(times, wind, min_windows=4):
    """
    Average the 10-minute wind means spanned by each PSD window.

    :type times: numpy.ndarray
    :param times: PSD window start times (Unix s)
    :type wind: dict
    :param wind: 10-minute wind means from hourly_wind
    :type min_windows: int
    :param min_windows: minimum number of 10-minute means needed (of 6)
    :rtype: numpy.ndarray
    :return: mean wind speed per PSD window, NaN if too little wind data
    """
    n = PSD_WINDOW_SEC // WIND_AVG_SEC
    out = np.full(len(times), np.nan)
    for i, t in enumerate(times):
        start = int(t) - int(t) % WIND_AVG_SEC
        vals = [wind[start + k * WIND_AVG_SEC] for k in range(n)
                if start + k * WIND_AVG_SEC in wind]
        if len(vals) >= min_windows:
            out[i] = np.mean(vals)
    return out


def bin_label(low, high):
    """
    Format a wind bin for the legend.

    :type low: float
    :param low: lower edge in m/s
    :type high: float
    :param high: upper edge in m/s, inf for an open bin
    :rtype: str
    :return: label, e.g. "1–2 m/s" or "> 5 m/s"
    """
    if np.isinf(high):
        return f"> {low:g} m/s"
    return f"{low:g}–{high:g} m/s"


def binned_spectra(station, comp, window, wind, edges, centers, args):
    """
    Median spectrum of one station and component per wind-speed bin.

    :type station: Station
    :param station: station to process
    :type comp: str
    :param comp: component code
    :type window: tuple of int
    :param window: (start, end) in Unix seconds
    :type wind: dict
    :param wind: 10-minute wind means from hourly_wind
    :type edges: list of float
    :param edges: wind bin edges in m/s; the last bin is open-ended
    :type centers: numpy.ndarray
    :param centers: frequency bin centers in Hz
    :type args: argparse.Namespace
    :param args: parsed command line arguments
    :rtype: tuple or None
    :return: (median dB [wind bins, freq bins], hours per wind bin), or
        None if the station lacks the component
    """
    data = channel_data(station, comp, window)
    if data is None or len(data[1]) == 0:
        return None
    binned = to_bins(data[0], data[2], centers, args.octave_fraction)
    speed = wind_for_segments(data[1], wind)
    uppers = list(edges[1:]) + [np.inf]
    medians = np.full((len(edges), centers.size), np.nan)
    counts = np.zeros(len(edges), dtype=int)
    for b, (low, high) in enumerate(zip(edges, uppers)):
        sel = (speed >= low) & (speed < high)
        counts[b] = sel.sum()
        if counts[b] >= args.min_hours:
            medians[b] = np.nanmedian(binned[sel], axis=0)
    return medians, counts


def plot(site, stations, wind, args):
    """
    Draw the wind-binned power-difference grid for one site and save it.

    :type site: str
    :param site: site code
    :type stations: list of Station
    :param stations: seismic stations at the site
    :type wind: dict
    :param wind: 10-minute wind means from hourly_wind
    :type args: argparse.Namespace
    :param args: parsed command line arguments
    """
    base = args.fontsize
    small = base - 2
    plt.rcParams.update({
        "font.size": base, "axes.labelsize": base,
        "axes.titlesize": base + 2, "xtick.labelsize": base - 1,
        "ytick.labelsize": base - 1, "legend.fontsize": small,
    })
    edges = args.bins
    centers = log_bins(args.fmin, args.fmax, args.octave_fraction)
    window = common_window(stations)
    rows = sorted(stations, key=lambda s: KEY_ORDER.index(station_key(s))
                  if station_key(s) in KEY_ORDER else 99)
    comps = args.components
    nrow, ncol = len(rows), len(comps)
    fig, axes = plt.subplots(
        nrow, ncol, sharex="col", sharey=True, squeeze=False,
        figsize=(args.panel_width * ncol, args.panel_height * nrow + 0.9),
        layout="constrained",
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=args.gap,
                                hspace=args.gap)

    ramp = plt.get_cmap(WIND_CMAP)
    colors = [CALM_COLOR] + [ramp(x) for x in
                             np.linspace(0.3, 1.0, len(edges) - 1)]
    uppers = list(edges[1:]) + [np.inf]
    handles = {}
    hour_counts = np.zeros(len(edges), dtype=int)
    for r, station in enumerate(rows):
        for c, comp in enumerate(comps):
            ax = axes[r, c]
            ax.axhline(0, color=CALM_COLOR, lw=1.0, zorder=1)
            result = binned_spectra(station, comp, window, wind, edges,
                                    centers, args)
            if result is None or np.isnan(result[0][0]).all():
                ax.text(0.5, 0.5, "no data" if result is None
                        else "too few calm hours", transform=ax.transAxes,
                        ha="center", va="center", fontsize=small,
                        color=MUTED_TEXT)
            else:
                medians, counts = result
                if comp == comps[0]:
                    hour_counts = np.maximum(hour_counts, counts)
                calm = medians[0]
                for b in range(1, len(edges)):
                    if np.isnan(medians[b]).all():
                        continue
                    h = ax.plot(centers, medians[b] - calm,
                                color=colors[b], lw=1.6, zorder=2 + b)
                    handles.setdefault(b, h[0])
            code = (station.sta if station.net == "XX"
                    else f"{station.net}.{station.sta}")
            ax.text(0.03, 0.95, code, transform=ax.transAxes, va="top",
                    fontsize=small, color=MUTED_TEXT)
            if r == 0:
                ax.set_title(comp)
            if c == 0:
                ax.set_ylabel(key_label(station_key(station)),
                              fontweight="bold")
            if r == nrow - 1:
                ax.set_xlabel("Frequency (Hz)")
            style_axis(ax, args)
    handles.setdefault(0, axes[0, 0].lines[0])

    # Legend: wind bins with the number of PSD windows in each (from the
    # station with the most), calm reference first
    order = sorted(handles)
    labels = {}
    for b in order:
        ref = "ref., " if b == 0 else ""
        labels[b] = (f"{bin_label(edges[b], uppers[b])} "
                     f"({ref}n={hour_counts[b]})")
    # Matplotlib fills legend columns top to bottom; reorder so entries
    # read left to right in increasing wind speed
    ncol = min(len(order), 3)
    nrows = int(np.ceil(len(order) / ncol))
    grid = [order[i:i + ncol] for i in range(0, len(order), ncol)]
    filled = [grid[r][c] for c in range(ncol) for r in range(nrows)
              if c < len(grid[r])]
    fig.legend([handles[b] for b in filled], [labels[b] for b in filled],
               loc="outside lower center", ncol=ncol, frameon=False,
               title=f"Wind speed at {site}{args.suffix} "
                     "(n = number of 1-h PSD windows)")

    start, end = (str(np.datetime64(int(t), "s"))[:10] for t in window)
    fig.suptitle(f"{site}: {start} to {end}", fontsize=base + 2)
    fig.supylabel(f"Power difference from {bin_label(edges[0], edges[1])}"
                  " (dB)", fontsize=base)

    outfile = Path(args.outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, dpi=args.dpi)
    plt.close(fig)
    print(f"Saved {outfile}")


def parse_args():
    """
    Parse command line arguments.

    :rtype: argparse.Namespace
    :return: parsed arguments
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--site", default="MITD",
                        help="Site to plot. Default: MITD")
    parser.add_argument("--spectra", default="SPECTRA",
                        help="Directory with *_all.npz files")
    parser.add_argument("--data", default="DATA",
                        help="MiniSEED directory with wind (LWS/LWD) files")
    parser.add_argument("--suffix", default="A",
                        help="Suffix of the anemometer station. Default: A")
    parser.add_argument("--bins", default="0,1,2,3,4,5",
                        help="Wind bin edges in m/s; the first bin is the "
                             "reference and the last is open-ended. "
                             "Default: 0,1,2,3,4,5")
    parser.add_argument("--min-hours", type=int, default=10,
                        help="Minimum PSD windows a wind bin needs to be "
                             "drawn. Default: 10")
    parser.add_argument("--components", default="Z,N,E",
                        help="Comma list of components. Default: Z,N,E")
    parser.add_argument("--exclude", default="",
                        help="Comma list of channels to leave out, as "
                             "NET.STA.LOC.CHA")
    parser.add_argument("--fmin", type=float, default=0.1,
                        help="Lowest frequency in Hz. Default: 0.1")
    parser.add_argument("--fmax", type=float, default=10.0,
                        help="Highest frequency in Hz. Default: 10")
    parser.add_argument("--octave-fraction", type=int, default=12,
                        help="Smoothing bins per octave. Default: 12")
    parser.add_argument("--panel-width", type=float, default=3.4,
                        help="Panel width in inches. Default: 3.4")
    parser.add_argument("--panel-height", type=float, default=2.0,
                        help="Panel height in inches. Default: 2.0")
    parser.add_argument("--gap", type=float, default=0.01,
                        help="Space between panels as a fraction of the "
                             "panel size. Default: 0.01")
    parser.add_argument("--fontsize", type=float, default=14,
                        help="Base font size in points. Default: 14")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--outfile", default=None,
                        help="Output image. Default: "
                             "FIGURES/wind_binned_<site>.png")
    args = parser.parse_args()

    args.site = args.site.upper()
    args.bins = [float(v) for v in args.bins.split(",")]
    if len(args.bins) < 2 or any(np.diff(args.bins) <= 0):
        parser.error("--bins needs at least two increasing edges")
    args.components = [c.strip().upper() for c in args.components.split(",")]
    args.exclude = [c for c in args.exclude.split(",") if c]
    if args.outfile is None:
        args.outfile = f"FIGURES/wind_binned_{args.site}.png"
    return args


def main():
    """Load spectra and wind for one site and draw the figure"""
    args = parse_args()
    by_site = load_stations(args.spectra, {args.site}, args.exclude)
    if args.site not in by_site:
        raise SystemExit(f"No NPZ files found for site {args.site}")
    stations = by_site[args.site]
    print(f"{args.site}: {', '.join(f'{s.net}.{s.sta}' for s in stations)}")
    wind = hourly_wind(args.site, args.data, args.suffix)
    plot(args.site, stations, wind, args)


if __name__ == "__main__":
    main()
