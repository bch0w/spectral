#!/usr/bin/env python
"""
Composite PPSD figure for one component (vertical by default) at every
station.

One row per site and one column per station type (Permanent, Covered,
Uncovered, Windscreen). Each panel is a probabilistic PSD: the hourly
acceleration PSDs inside the site's common time window are smoothed onto
1/12-octave frequency bins (as in plot_permanent_relative.py), and for each
frequency the share of hours falling in each 1 dB power bin is shown in
color. The Peterson NLNM/NHNM are drawn on top. The colormap is set with
--cmap (any Matplotlib colormap name, or the default "blue").

Sites are ordered by the start of their common time window, then by name.
All sites are included by default, KAED and RND too.

Example::

    python plot_ppsd_composite.py
    python plot_ppsd_composite.py --component N \
        --outfile FIGURES/ppsd_composite_N.png
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

from plot_permanent_relative import (  # noqa: E402
    KEY_ORDER, MUTED_TEXT, channel_data, common_window, key_label,
    load_stations, log_bins, station_key, style_axis, to_bins,
)


# Default sequential single-hue ramp (reference palette blue, light ->
# dark), selected with --cmap blue
BLUE_RAMP = LinearSegmentedColormap.from_list("blue", [
    "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b",
])
MODEL_COLOR = "r"
# Figure height (inches) added for the horizontal colorbar at the bottom
COLORBAR_HEIGHT = 0.8


def get_cmap(name):
    """
    Return the colormap to use for probabilities, with zero probability
    left transparent so empty bins read as background.

    :type name: str
    :param name: "blue" for the default ramp, or any Matplotlib colormap
        name such as viridis, cividis or magma_r
    :rtype: matplotlib.colors.Colormap
    :return: colormap copy with its "under" color set to transparent
    """
    if name == "blue":
        cmap = BLUE_RAMP.copy()
    else:
        try:
            cmap = matplotlib.colormaps[name].copy()
        except KeyError:
            raise SystemExit(f"Unknown colormap: {name}")
    cmap.set_under("none")
    return cmap


def model_in_freq(model):
    """
    Convert a Peterson model from period to frequency, sorted by frequency.

    :type model: tuple of numpy.ndarray
    :param model: (periods in s, dB)
    :rtype: tuple of numpy.ndarray
    :return: (frequencies in Hz, dB)
    """
    periods, db = (np.asarray(a, dtype=float) for a in model)
    freqs = 1.0 / periods
    order = np.argsort(freqs)
    return freqs[order], db[order]


def ppsd_panel(station, comp, window, centers, db_edges, args):
    """
    Build the probability density for one station and component.

    :type station: Station
    :param station: station to plot
    :type comp: str
    :param comp: component code, e.g. Z
    :type window: tuple of int
    :param window: (start, end) in Unix seconds
    :type centers: numpy.ndarray
    :param centers: frequency bin centers in Hz
    :type db_edges: numpy.ndarray
    :param db_edges: power bin edges in dB
    :type args: argparse.Namespace
    :param args: parsed command line arguments
    :rtype: numpy.ndarray or None
    :return: probability, shape [dB bins, freq bins], or None if the
        station lacks the component
    """
    data = channel_data(station, comp, window)
    if data is None or len(data[1]) == 0:
        return None
    binned = to_bins(data[0], data[2], centers, args.octave_fraction)
    prob = np.full((db_edges.size - 1, centers.size), np.nan)
    for i in range(centers.size):
        column = binned[:, i]
        column = column[np.isfinite(column)]
        if column.size:
            counts, _ = np.histogram(column, bins=db_edges)
            prob[:, i] = counts / column.size
    return prob


def edges_from_centers(centers, octave_fraction):
    """
    Return log-spaced bin edges around bin centers.

    :type centers: numpy.ndarray
    :param centers: bin centers in Hz
    :type octave_fraction: int
    :param octave_fraction: bins per octave
    :rtype: numpy.ndarray
    :return: bin edges in Hz, one more than the number of centers
    """
    half = 2 ** (0.5 / octave_fraction)
    return np.append(centers / half, centers[-1] * half)


def plot(by_site, args):
    """
    Draw the composite PPSD grid and save it.

    :type by_site: dict
    :param by_site: site code -> list of Station, in row order
    :type args: argparse.Namespace
    :param args: parsed command line arguments
    """
    base = args.fontsize
    small = base - 2
    plt.rcParams.update({
        "font.size": base, "axes.labelsize": base,
        "axes.titlesize": base + 2, "xtick.labelsize": base - 1,
        "ytick.labelsize": base - 1,
    })
    columns = [k for k in KEY_ORDER
               if any(station_key(s) == k for st in by_site.values()
                      for s in st)]
    sites = list(by_site)
    nrow, ncol = len(sites), len(columns)
    # Extra height holds the horizontal colorbar under the panels, so the
    # panels keep their height; with no colorbar at the side they use the
    # full figure width
    fig, axes = plt.subplots(
        nrow, ncol, sharex="col", sharey=True, squeeze=False,
        figsize=(args.panel_width * ncol,
                 args.panel_height * nrow + COLORBAR_HEIGHT),
        layout="constrained",
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=args.gap,
                                hspace=args.gap)

    centers = log_bins(args.fmin, args.fmax, args.octave_fraction)
    f_edges = edges_from_centers(centers, args.octave_fraction)
    db_edges = np.arange(args.ylim[0], args.ylim[1] + args.db_step,
                         args.db_step)
    any_ppsd = next(iter(by_site[sites[0]][0].ppsds.values()))
    nlnm = model_in_freq(any_ppsd.get_nlnm())
    nhnm = model_in_freq(any_ppsd.get_nhnm())

    cmap = get_cmap(args.cmap)
    image = None
    for row, site in enumerate(sites):
        stations = {station_key(s): s for s in by_site[site]}
        window = common_window(by_site[site])
        for col, key in enumerate(columns):
            ax = axes[row, col]
            station = stations.get(key)
            panel = (ppsd_panel(station, args.component, window, centers,
                                db_edges, args)
                     if station is not None else None)
            if panel is None:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=small,
                        color=MUTED_TEXT)
            else:
                image = ax.pcolormesh(f_edges, db_edges, panel,
                                      cmap=cmap, vmin=1e-9,
                                      vmax=args.pmax, shading="flat",
                                      rasterized=True)
                # White backing keeps the code readable over the NLNM line;
                # XX stations are labeled by station code only
                code = (station.sta if station.net == "XX"
                        else f"{station.net}.{station.sta}")
                ax.text(0.05, 0.05, code,
                        transform=ax.transAxes, ha="left", va="bottom",
                        fontsize=small, color=MUTED_TEXT,
                        bbox=dict(facecolor="white", edgecolor="none",
                                  alpha=0.85, pad=1.5))
            for model in (nlnm, nhnm):
                ax.plot(*model, color=MODEL_COLOR, lw=1.0, ls="--")

            if row == 0:
                ax.set_title(key_label(key))
            if col == 0:
                ax.set_ylabel(site, fontweight="bold")
                # Short date form so the label fits inside the panel
                start, end = (str(np.datetime64(int(t), "s"))[:10]
                              for t in window)
                end = end[5:] if end[:4] == start[:4] else end
                ax.text(0.03, 0.95, f"{start} – {end}",
                        transform=ax.transAxes, va="top", fontsize=small,
                        color=MUTED_TEXT,
                        bbox=dict(facecolor="white", edgecolor="none",
                                  alpha=0.85, pad=1.5))
            if row == nrow - 1:
                ax.set_xlabel("Frequency (Hz)")
            style_axis(ax, args)
            ax.set_ylim(*args.ylim)
            # No grid over the density (minor lines show as seams)
            ax.grid(False, which="both")

    fig.supylabel(f"{args.component}-component acceleration PSD "
                  r"(dB rel. 1 (m/s$^2$)$^2$/Hz)", fontsize=base)
    if image is not None:
        bar = fig.colorbar(image, ax=axes, location="bottom", shrink=0.5,
                           aspect=30, pad=0.01, extend="max")
        # A few short ticks so the horizontal labels don't run together
        bar.set_ticks(np.linspace(0, args.pmax, 4))
        bar.ax.xaxis.set_major_formatter("{x:.2f}")
        bar.set_label("Probability")
        bar.outline.set_visible(False)

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
    parser.add_argument("--component", default="Z",
                        help="Component to plot. Default: Z")
    parser.add_argument("--sites", default=None,
                        help="Comma list of sites. Default: all")
    parser.add_argument("--exclude-sites", default="",
                        help="Comma list of sites to leave out. "
                             "Default: none")
    parser.add_argument("--exclude", default="",
                        help="Comma list of channels to leave out, as "
                             "NET.STA.LOC.CHA")
    parser.add_argument("--fmin", type=float, default=0.1,
                        help="Lowest frequency in Hz. Default: 0.1")
    parser.add_argument("--fmax", type=float, default=10.0,
                        help="Highest frequency in Hz. Default: 10")
    parser.add_argument("--octave-fraction", type=int, default=12,
                        help="Frequency bins per octave. Default: 12")
    parser.add_argument("--ylim", type=float, nargs=2,
                        default=(-180.0, -80.0),
                        help="Power range in dB. Default: -180 -80")
    parser.add_argument("--db-step", type=float, default=1.0,
                        help="Power bin width in dB. Default: 1")
    parser.add_argument("--cmap", default="blue",
                        help="Colormap for probability: 'blue' (default "
                             "single-hue ramp) or any Matplotlib colormap "
                             "name, e.g. viridis, cividis, magma_r")
    parser.add_argument("--pmax", type=float, default=0.15,
                        help="Probability at the top of the color scale. "
                             "Default: 0.15")
    # Wider panels make a wider figure that is shrunk more to fit the
    # page, so the default font is larger to keep printed text ~9 pt
    parser.add_argument("--panel-width", type=float, default=3.0,
                        help="Panel width in inches. Default: 3.0")
    parser.add_argument("--panel-height", type=float, default=1.5,
                        help="Panel height in inches. Default: 1.5")
    parser.add_argument("--gap", type=float, default=0.01,
                        help="Space between panels as a fraction of the "
                             "panel size. Default: 0.01")
    parser.add_argument("--fontsize", type=float, default=17,
                        help="Base font size in points. Default: 17")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--outfile", default=None,
                        help="Output image. Default: "
                             "FIGURES/ppsd_composite_<component>.png")
    args = parser.parse_args()

    args.component = args.component.upper()
    args.sites = set(args.sites.split(",")) if args.sites else None
    args.exclude_sites = {s for s in args.exclude_sites.split(",") if s}
    args.exclude = [c for c in args.exclude.split(",") if c]
    if args.outfile is None:
        args.outfile = f"FIGURES/ppsd_composite_{args.component}.png"
    return args


def main():
    """Load the spectra and draw the composite figure"""
    args = parse_args()
    by_site = load_stations(args.spectra, args.sites, args.exclude)
    for site in args.exclude_sites & set(by_site):
        del by_site[site]
    if not by_site:
        raise SystemExit("No matching NPZ files found")

    # Same row order as the comparison figures: start day, then name
    order = sorted(by_site,
                   key=lambda s: (common_window(by_site[s])[0] // 86400, s))
    by_site = {site: by_site[site] for site in order}
    for site, stations in by_site.items():
        names = ", ".join(f"{s.net}.{s.sta}" for s in stations)
        print(f"{site}: {names}")
    plot(by_site, args)


if __name__ == "__main__":
    main()
