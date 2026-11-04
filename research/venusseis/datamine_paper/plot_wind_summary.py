#!/usr/bin/env python
"""
Summary of wind-induced noise at every site, one curve per installation.

Rows are sites and columns are components (Z, N, E). Each panel has one
curve per station (Permanent, Covered, Uncovered, Windscreen): the
station's median spectrum during windy hours minus its own median spectrum
during calm hours, using the wind measured at the site's anemometer
(<SITE>A). Positive values mean more noise when it is windy.

Windy and calm hours use the same thresholds at every site (default: calm
< 1 m/s, windy >= 4 m/s), so curves compare directly. Sites without enough
calm or windy hours are left out, with a note printed. PSDs are smoothed
and limited to each site's common window as in plot_permanent_relative.py.
Wind is read with plot_wind.py (faulty samples removed), averaged in
10-minute windows, and each PSD hour gets the mean of the 10-minute windows
it spans.

Example::

    python plot_wind_summary.py
    python plot_wind_summary.py --windy 5 --sites K24K,MITD
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from plot_permanent_relative import (  # noqa: E402
    KEY_ORDER, MUTED_TEXT, PERMANENT_COLOR, SUFFIX_COLORS, channel_data,
    common_window, key_label, load_stations, log_bins, station_key,
    style_axis, to_bins,
)
from plot_wind import find_wind_files, load_station  # noqa: E402


# Wind is averaged in 10-minute windows; each PSD covers one hour
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
    :return: (median dB [wind bins, freq bins], PSD windows per wind bin),
        or None if the station lacks the component; bins with fewer than
        args.min_hours windows are NaN
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


def site_curves(stations, wind, comps, centers, args):
    """
    Windy-minus-calm spectra for every station and component at a site.

    :type stations: list of Station
    :param stations: stations at the site
    :type wind: dict
    :param wind: 10-minute wind means from hourly_wind
    :type comps: list of str
    :param comps: component codes
    :type centers: numpy.ndarray
    :param centers: frequency bin centers in Hz
    :type args: argparse.Namespace
    :param args: parsed command line arguments
    :rtype: tuple
    :return: (dict (key, comp) -> dB curve, calm windows, windy windows,
        common time window); window counts are the largest over stations
    """
    # Bins: [0, calm), [calm, windy), [windy, inf); the middle is unused
    edges = [0.0, args.calm, args.windy]
    window = common_window(stations)
    curves, n_calm, n_windy = {}, 0, 0
    for station in stations:
        for comp in comps:
            result = binned_spectra(station, comp, window, wind, edges,
                                    centers, args)
            if result is None:
                continue
            medians, counts = result
            n_calm = max(n_calm, counts[0])
            n_windy = max(n_windy, counts[2])
            if np.isnan(medians[0]).all() or np.isnan(medians[2]).all():
                continue
            curves[(station_key(station), comp)] = medians[2] - medians[0]
    return curves, n_calm, n_windy, window


def plot(results, args):
    """
    Draw the summary grid and save it.

    :type results: dict
    :param results: site -> output of site_curves, in row order
    :type args: argparse.Namespace
    :param args: parsed command line arguments
    """
    base = args.fontsize
    small = base - 2
    plt.rcParams.update({
        "font.size": base, "axes.labelsize": base,
        "axes.titlesize": base + 2, "xtick.labelsize": base - 1,
        "ytick.labelsize": base - 1, "legend.fontsize": base,
    })
    comps = args.components
    sites = list(results)
    nrow, ncol = len(sites), len(comps)
    fig, axes = plt.subplots(
        nrow, ncol, sharex="col", sharey="row", squeeze=False,
        figsize=(args.panel_width * ncol, args.panel_height * nrow),
        layout="constrained",
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=args.gap,
                                hspace=args.gap)
    centers = log_bins(args.fmin, args.fmax, args.octave_fraction)
    colors = {"Permanent": PERMANENT_COLOR, **SUFFIX_COLORS}

    handles = {}
    for r, site in enumerate(sites):
        curves, n_calm, n_windy, window = results[site]
        lo, hi = 0.0, 0.0
        for c, comp in enumerate(comps):
            ax = axes[r, c]
            ax.axhline(0, color="#8a8984", lw=1.0, zorder=1)
            for key in KEY_ORDER:
                curve = curves.get((key, comp))
                if curve is None:
                    continue
                h = ax.plot(centers, curve, color=colors[key], lw=1.6,
                            zorder=3 if key == "Permanent" else 2)
                handles.setdefault(key, h[0])
                lo = min(lo, np.nanmin(curve))
                hi = max(hi, np.nanmax(curve))
            if r == 0:
                ax.set_title(comp)
            if c == 0:
                ax.set_ylabel(site, fontweight="bold")
                start, end = (str(np.datetime64(int(t), "s"))[:10]
                              for t in window)
                end = end[5:] if end[:4] == start[:4] else end
                ax.text(0.03, 0.95, f"{start} – {end}\n"
                        f"n = {n_calm} calm, {n_windy} windy",
                        transform=ax.transAxes, va="top", fontsize=small,
                        color=MUTED_TEXT)
            if r == nrow - 1:
                ax.set_xlabel("Frequency (Hz)")
            style_axis(ax, args)
        # Row y-limits from the data, with headroom for the label
        span = (hi - lo) or 1.0
        axes[r, 0].set_ylim(lo - 0.05 * span, hi + 0.35 * span)

    order = [k for k in KEY_ORDER if k in handles]
    fig.legend([handles[k] for k in order], [key_label(k) for k in order],
               loc="outside upper center", ncol=len(order), frameon=False)
    fig.supylabel(f"Power difference, wind ≥ {args.windy:g} m/s minus "
                  f"< {args.calm:g} m/s (dB)", fontsize=base)

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
    parser.add_argument("--spectra", default="SPECTRA",
                        help="Directory with *_all.npz files")
    parser.add_argument("--data", default="DATA",
                        help="MiniSEED directory with wind (LWS/LWD) files")
    parser.add_argument("--suffix", default="A",
                        help="Suffix of the anemometer station. Default: A")
    parser.add_argument("--sites", default=None,
                        help="Comma list of sites. Default: all with "
                             "enough calm and windy hours")
    parser.add_argument("--calm", type=float, default=1.0,
                        help="Calm hours have wind below this (m/s). "
                             "Default: 1")
    parser.add_argument("--windy", type=float, default=4.0,
                        help="Windy hours have wind at or above this "
                             "(m/s). Default: 4")
    parser.add_argument("--min-hours", type=int, default=20,
                        help="Minimum PSD windows needed in both the calm "
                             "and windy groups. Default: 20")
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
    parser.add_argument("--outfile", default="FIGURES/wind_summary.png",
                        help="Output image. Default: "
                             "FIGURES/wind_summary.png")
    args = parser.parse_args()

    if not 0 < args.calm <= args.windy:
        parser.error("need 0 < --calm <= --windy")
    args.sites = set(args.sites.split(",")) if args.sites else None
    args.components = [c.strip().upper() for c in args.components.split(",")]
    args.exclude = [c for c in args.exclude.split(",") if c]
    return args


def main():
    """Compute windy-minus-calm spectra for every site and plot them"""
    args = parse_args()
    by_site = load_stations(args.spectra, args.sites, args.exclude)
    if not by_site:
        raise SystemExit("No matching NPZ files found")

    # Same row order as the other figures: start day, then name
    order = sorted(by_site,
                   key=lambda s: (common_window(by_site[s])[0] // 86400, s))
    centers = log_bins(args.fmin, args.fmax, args.octave_fraction)
    results = {}
    for site in order:
        try:
            wind = hourly_wind(site, args.data, args.suffix)
        except SystemExit as err:
            print(f"{site}: skipped ({err})")
            continue
        out = site_curves(by_site[site], wind, args.components, centers,
                          args)
        curves, n_calm, n_windy, _ = out
        if not curves:
            print(f"{site}: skipped, too few hours "
                  f"(calm {n_calm}, windy {n_windy}; need "
                  f"{args.min_hours} each)")
            continue
        print(f"{site}: {n_calm} calm and {n_windy} windy PSD windows")
        results[site] = out
    if not results:
        raise SystemExit("No site has enough calm and windy hours")
    plot(results, args)


if __name__ == "__main__":
    main()
