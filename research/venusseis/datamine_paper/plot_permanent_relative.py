#!/usr/bin/env python
"""
Plot noise at co-located XX stations relative to the site's permanent
station.

One row per site and one column per component. Each XX station's PSD is
compared with the permanent station hour by hour: both are smoothed onto the
same log-spaced frequency bins, hours that both stations recorded are
paired, and the difference (XX minus permanent, in dB) is taken for each
pair. Each panel shows the median of those hourly differences with a shaded
interquartile band; 0 dB means the same noise as the permanent station.
All stations at a site are limited to their common time window.

All sites go in one figure, ordered by the start of their common time window
and then by site name. Panels in a row share y-limits; rows can differ.
KAED and RND are left out by default (see --exclude-sites);
they are handled in a separate figure (plot_covered_relative.py), which
reuses this code with --reference C.

Reads the combined ``NET_STA_CHA_all.npz`` files written by
``run_psd_viz.py`` / ``psd_tools.py compute``.

Example::

    python plot_permanent_relative.py
    python plot_permanent_relative.py --sites MITD --components Z,H
"""
import argparse
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from watc_ppsd import WelchPPSD  # noqa: E402


# Network holding the co-located sub-stations of each permanent station
SUB_NETWORK = "XX"

# NPZ names written by psd_tools.py, e.g. XX_MITDC_HHZ_all.npz
NPZ_RE = re.compile(
    r"^(?P<net>[A-Z0-9]+)_(?P<sta>[A-Z0-9]+)_(?:(?P<loc>[A-Z0-9]+)_)?"
    r"(?P<cha>[A-Z0-9]{3})_all\.npz$"
)

# Light-mode categorical palette (validated default, first three slots,
# which pass all-pairs CVD checks). The permanent station is the 0 dB line.
SUFFIX_COLORS = {"C": "#2a78d6", "U": "#eb6834", "W": "#1baf7a"}
SUFFIX_LABELS = {"C": "Covered", "U": "Uncovered", "W": "Windscreen"}

# Sites handled in their own figure: KAED has a broken permanent reference,
# RND shows a shared high-frequency offset that points at the reference
DEFAULT_EXCLUDED_SITES = "KAED,RND"
FALLBACK_COLORS = ["#eda100", "#e87ba4", "#008300", "#6250d6"]
PERMANENT_COLOR = "#3d3d3a"
MUTED_TEXT = "#52514e"

# Legend order, by station key ("Permanent" or an XX suffix)
KEY_ORDER = ["Permanent", "C", "U", "W"]


def station_key(station):
    """
    Return the key used for a station's color, label and legend entry.

    :type station: Station
    :param station: station to name
    :rtype: str
    :return: "Permanent" or the XX suffix, e.g. "C"
    """
    return "Permanent" if station.is_permanent else station.suffix


def key_label(key):
    """
    Return the legend label for a station key.

    :type key: str
    :param key: "Permanent" or an XX suffix
    :rtype: str
    :return: label, e.g. "Covered"
    """
    return SUFFIX_LABELS.get(key, key)


class Station:
    """
    PSD data for one station at one site, keyed by component.

    :type net: str
    :param net: network code
    :type sta: str
    :param sta: station code
    :type site: str
    :param site: code of the site this station belongs to, e.g. MITD
    """
    def __init__(self, net, sta, site):
        self.net = net
        self.sta = sta
        self.site = site
        self.ppsds = {}

    @property
    def is_permanent(self):
        """True for the site's permanent (non-XX) station"""
        return self.net != SUB_NETWORK

    @property
    def suffix(self):
        """Suffix that distinguishes an XX station, e.g. C for MITDC"""
        return "" if self.is_permanent else self.sta[len(self.site):]


def site_of(net, sta, permanent_names):
    """
    Map a station to its site code. XX stations are named after the
    permanent station plus a suffix, e.g. MITDC, MITDU, MITDW -> MITD.

    :type net: str
    :param net: network code
    :type sta: str
    :param sta: station code
    :type permanent_names: set of str
    :param permanent_names: codes of all permanent (non-XX) stations
    :rtype: str
    :return: site code
    """
    if net != SUB_NETWORK:
        return sta
    for name in sorted(permanent_names, key=len, reverse=True):
        if sta.startswith(name) and sta != name:
            return name
    return sta[:-1]


def load_stations(spectra_dir, sites=None, exclude=()):
    """
    Load every combined NPZ in a directory and group the stations by site.

    :type spectra_dir: str
    :param spectra_dir: directory with ``*_all.npz`` files
    :type sites: set of str
    :param sites: optional site codes to keep; None keeps all
    :type exclude: iterable of str
    :param exclude: channel codes to skip, as NET.STA.LOC.CHA
    :rtype: dict
    :return: site code -> list of Station
    """
    exclude = set(exclude)
    found = []
    for path in sorted(Path(spectra_dir).glob("*_all.npz")):
        match = NPZ_RE.match(path.name)
        if not match:
            continue
        sncl = (f"{match['net']}.{match['sta']}.{match['loc'] or ''}."
                f"{match['cha']}")
        if sncl in exclude:
            print(f"Excluding {sncl}")
            continue
        found.append((match["net"], match["sta"], match["cha"], path))

    permanent = {sta for net, sta, _, _ in found if net != SUB_NETWORK}
    stations = {}
    for net, sta, cha, path in found:
        site = site_of(net, sta, permanent)
        if sites and site not in sites:
            continue
        key = (net, sta)
        if key not in stations:
            stations[key] = Station(net, sta, site)
        stations[key].ppsds[cha[-1]] = WelchPPSD.from_npz(path)

    by_site = defaultdict(list)
    for station in stations.values():
        by_site[station.site].append(station)
    for members in by_site.values():
        # Permanent station first, then XX stations by suffix
        members.sort(key=lambda s: (not s.is_permanent, s.suffix))
    return dict(by_site)


def segment_times(ppsd):
    """
    Return the start times of a PPSD's segments as whole Unix seconds.

    :type ppsd: watc_ppsd.WelchPPSD
    :param ppsd: PPSD to read
    :rtype: numpy.ndarray
    :return: segment start times
    """
    return np.array([round(t.timestamp) for t in ppsd.current_times_used],
                    dtype=np.int64)


def common_window(stations):
    """
    Find the time span covered by every channel of every station at a site.

    :type stations: list of Station
    :param stations: stations at one site
    :rtype: tuple of int
    :return: (start, end) in Unix seconds
    """
    spans = [segment_times(p) for s in stations for p in s.ppsds.values()]
    start = max(t.min() for t in spans)
    end = min(t.max() for t in spans)
    if start >= end:
        raise ValueError(f"No common time window at {stations[0].site}")
    return start, end


def log_bins(fmin, fmax, octave_fraction):
    """
    Build log-spaced frequency bin centers.

    :type fmin: float
    :param fmin: lowest bin center in Hz
    :type fmax: float
    :param fmax: highest bin center in Hz
    :type octave_fraction: int
    :param octave_fraction: bins per octave, e.g. 12 for 1/12 octave
    :rtype: numpy.ndarray
    :return: bin centers in Hz
    """
    n_bins = int(np.ceil(np.log2(fmax / fmin) * octave_fraction)) + 1
    return fmin * 2 ** (np.arange(n_bins) / octave_fraction)


def to_bins(freqs, values, centers, octave_fraction):
    """
    Smooth PSD segments onto log-spaced frequency bins by averaging the dB
    values in each bin. Bins narrower than the frequency spacing are
    interpolated; bins above the highest frequency are NaN.

    :type freqs: numpy.ndarray
    :param freqs: frequencies in Hz
    :type values: numpy.ndarray
    :param values: dB values, shape [segments, frequencies]
    :type centers: numpy.ndarray
    :param centers: bin centers in Hz
    :type octave_fraction: int
    :param octave_fraction: bins per octave used to build ``centers``
    :rtype: numpy.ndarray
    :return: binned dB values, shape [segments, bins]
    """
    half = 2 ** (0.5 / octave_fraction)
    out = np.full((values.shape[0], centers.size), np.nan)
    valid = freqs > 0
    log_f = np.log10(freqs[valid])
    order = np.argsort(log_f)
    log_f = log_f[order]
    vv = values[:, valid][:, order]
    for i, fc in enumerate(centers):
        if fc > freqs.max():
            continue
        in_bin = (freqs >= fc / half) & (freqs < fc * half)
        if in_bin.any():
            out[:, i] = values[:, in_bin].mean(axis=1)
        else:
            # Bin narrower than the frequency spacing: interpolate
            j = int(np.clip(np.searchsorted(log_f, np.log10(fc)), 1,
                            log_f.size - 1))
            w = (np.log10(fc) - log_f[j - 1]) / (log_f[j] - log_f[j - 1])
            out[:, i] = (1 - w) * vv[:, j - 1] + w * vv[:, j]
    return out


def channel_data(station, comp, window):
    """
    Return one station's PSD segments for a component inside a time window.
    Component H combines N and E by averaging their power (not their dB) in
    each hour both channels recorded.

    :type station: Station
    :param station: station to read
    :type comp: str
    :param comp: Z, N, E or H
    :type window: tuple of int
    :param window: (start, end) in Unix seconds
    :rtype: tuple or None
    :return: (frequencies, segment times, dB array [segments, frequencies]),
        or None if the station lacks the component
    """
    def one(ppsd):
        times = segment_times(ppsd)
        keep = (times >= window[0]) & (times <= window[1])
        return (np.asarray(ppsd.psd_frequencies, dtype=float), times[keep],
                np.asarray(ppsd.psd_values)[keep])

    if comp != "H":
        return one(station.ppsds[comp]) if comp in station.ppsds else None
    if not {"N", "E"} <= set(station.ppsds):
        return None
    f_n, t_n, v_n = one(station.ppsds["N"])
    f_e, t_e, v_e = one(station.ppsds["E"])
    if not np.allclose(f_n, f_e):
        raise ValueError(f"N and E grids differ for {station.sta}")
    shared, i_n, i_e = np.intersect1d(t_n, t_e, return_indices=True)
    power = (10 ** (v_n[i_n] / 10) + 10 ** (v_e[i_e] / 10)) / 2
    return f_n, shared, 10 * np.log10(power)


def paired_difference(data, ref, centers, octave_fraction):
    """
    Compute hourly differences between a station and the permanent station
    on common frequency bins, using only hours both recorded.

    :type data: tuple
    :param data: XX station output of channel_data
    :type ref: tuple
    :param ref: permanent station output of channel_data
    :type centers: numpy.ndarray
    :param centers: frequency bin centers in Hz
    :type octave_fraction: int
    :param octave_fraction: bins per octave used to build ``centers``
    :rtype: numpy.ndarray or None
    :return: dB differences, shape [paired hours, bins], or None if the
        stations share no hours
    """
    shared, i_sta, i_ref = np.intersect1d(data[1], ref[1],
                                          return_indices=True)
    if shared.size == 0:
        return None
    sta_binned = to_bins(data[0], data[2][i_sta], centers, octave_fraction)
    ref_binned = to_bins(ref[0], ref[2][i_ref], centers, octave_fraction)
    return sta_binned - ref_binned


def style_axis(ax, args, last_col):
    """
    Apply the shared axis styling. Tick labels are plain numbers (0.1, 1,
    10); only the last column labels the right edge so that labels of
    neighboring panels do not collide.

    :type ax: matplotlib.axes.Axes
    :param ax: axis to style
    :type args: argparse.Namespace
    :param args: parsed command line arguments
    :type last_col: bool
    :param last_col: True for panels in the right-most column
    """
    ax.set_xscale("log")
    ax.set_xlim(args.fmin, args.fmax)
    decades = 10.0 ** np.arange(np.ceil(np.log10(args.fmin)),
                                np.floor(np.log10(args.fmax)) + 1)
    labels = [f"{d:g}" for d in decades]
    if not last_col and np.isclose(decades[-1], args.fmax):
        labels[-1] = ""
    ax.set_xticks(decades, labels)
    ax.grid(True, which="both", lw=0.25, color="#d8d7d2")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def plot(by_site, components, outfile, args):
    """
    Draw the grid of permanent-relative spectra and save it.

    :type by_site: dict
    :param by_site: site code -> list of Station
    :type components: list of str
    :param components: component codes per column; H means combined N/E
    :type outfile: str
    :param outfile: output image path
    :type args: argparse.Namespace
    :param args: parsed command line arguments
    """
    sites = list(by_site)
    nrow, ncol = len(sites), len(components)
    base = args.fontsize
    small = base - 2
    plt.rcParams.update({
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base + 2,
        "xtick.labelsize": base - 1,
        "ytick.labelsize": base - 1,
        "legend.fontsize": base,
    })
    # Constrained layout packs panels with minimal padding; inner tick
    # labels are already hidden by the shared axes
    fig, axes = plt.subplots(
        nrow, ncol, sharex="col", sharey="row", squeeze=False,
        figsize=(args.panel_width * ncol, args.panel_height * nrow),
        layout="constrained",
    )
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.02, wspace=args.gap,
                                hspace=args.gap)
    centers = log_bins(args.fmin, args.fmax, args.octave_fraction)

    ref_key = args.reference
    ref_label = key_label(ref_key)
    ref_color = (PERMANENT_COLOR if ref_key == "Permanent"
                 else SUFFIX_COLORS.get(ref_key, PERMANENT_COLOR))

    handles = {}
    for row, site in enumerate(sites):
        data_range = [0.0, 0.0]
        stations = by_site[site]
        reference = next((s for s in stations if station_key(s) == ref_key),
                         None)
        window = common_window(stations)
        for col, comp in enumerate(components):
            ax = axes[row, col]
            h = ax.axhline(0, color=ref_color, lw=1.2, zorder=1)
            handles.setdefault(ref_key, h)

            ref = (channel_data(reference, comp, window)
                   if reference is not None else None)
            if ref is None:
                ax.text(0.5, 0.5, f"no {ref_label.lower()}-station "
                        "reference", transform=ax.transAxes, ha="center",
                        va="center", fontsize=small, color=MUTED_TEXT)

            spare = iter(FALLBACK_COLORS)
            for station in stations:
                if station is reference or ref is None:
                    continue
                data = channel_data(station, comp, window)
                if data is None:
                    continue
                diff = paired_difference(data, ref, centers,
                                         args.octave_fraction)
                if diff is None:
                    continue
                q25, q50, q75 = np.nanpercentile(diff, (25, 50, 75),
                                                 axis=0)
                key = station_key(station)
                if station.is_permanent:
                    color = PERMANENT_COLOR
                else:
                    color = SUFFIX_COLORS.get(key) or next(spare)
                ax.fill_between(centers, q25, q75, color=color, alpha=0.15,
                                lw=0, zorder=2)
                h = ax.plot(centers, q50, color=color, lw=1.6, zorder=3)
                handles.setdefault(key, h[0])
                data_range[0] = min(data_range[0], np.nanmin(q25))
                data_range[1] = max(data_range[1], np.nanmax(q75))

            if row == 0:
                ax.set_title("N+E" if comp == "H" else comp)
            if col == 0:
                start, end = (np.datetime64(int(t), "s") for t in window)
                # Site name per row; the shared "dB rel." label is added
                # once for the figure so short rows don't crowd
                ax.set_ylabel(site, fontweight="bold")
                ax.text(0.02, 0.97,
                        f"{str(start)[:10]} to {str(end)[:10]}",
                        transform=ax.transAxes, va="top", fontsize=small,
                        color=MUTED_TEXT)
            if row == nrow - 1:
                ax.set_xlabel("Frequency (Hz)")
            style_axis(ax, args, last_col=col == ncol - 1)

        # Panels in a row share y-limits fitted to that row's data (always
        # including 0 dB), with extra headroom for the date label
        if args.ylim:
            axes[row, 0].set_ylim(*args.ylim)
        else:
            span = (data_range[1] - data_range[0]) or 1.0
            axes[row, 0].set_ylim(data_range[0] - 0.05 * span,
                                  data_range[1] + 0.15 * span)

    order = [k for k in KEY_ORDER if k in handles]
    order += sorted(k for k in handles if k not in KEY_ORDER)
    labels = [key_label(k) for k in order]

    # One y-axis label for the whole figure, outside the site names
    fig.supylabel(f"dB rel. {ref_label.lower()} station", fontsize=base)

    # Legend in its own strip above the panels (constrained layout reserves
    # the space); wrap onto more rows when the figure is narrow
    entry_width = 1.2 * base / 10
    per_row = max(1, min(len(order),
                         int(fig.get_figwidth() / entry_width)))
    fig.legend([handles[k] for k in order], labels,
               loc="outside upper center", ncol=per_row, frameon=False,
               borderaxespad=0.2)

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, dpi=args.dpi)
    plt.close(fig)
    print(f"Saved {outfile}")


def parse_args(description=__doc__, **defaults):
    """
    Parse command line arguments.

    :type description: str
    :param description: help text shown by --help
    :type defaults: dict
    :param defaults: overrides for argument defaults, so other figure
        scripts can reuse this parser, e.g. ``reference="C"``
    :rtype: argparse.Namespace
    :return: parsed arguments
    """
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--reference", default="permanent",
                        help="Station each curve is compared with: "
                             "'permanent' or an XX suffix such as C "
                             "(Covered). Default: permanent")
    parser.add_argument("--omit", default="",
                        help="Comma list of stations to leave out at every "
                             "site: 'permanent' and/or XX suffixes, e.g. "
                             "permanent,W. Default: none")
    parser.add_argument("--spectra", default="SPECTRA",
                        help="Directory with *_all.npz files")
    parser.add_argument("--sites", default=None,
                        help="Comma list of sites, e.g. MITD,MDM. "
                             "Default: all")
    parser.add_argument("--components", default="Z,N,E",
                        help="Comma list of columns from Z, N, E and H "
                             "(H = power mean of N and E). Default: Z,N,E")
    parser.add_argument("--exclude", default="",
                        help="Comma list of channels to leave out, as "
                             "NET.STA.LOC.CHA, e.g. HV.KAED..EHZ")
    parser.add_argument("--fmin", type=float, default=0.1,
                        help="Lowest frequency in Hz. Default: 0.1")
    parser.add_argument("--fmax", type=float, default=10.0,
                        help="Highest frequency in Hz. Default: 10")
    parser.add_argument("--octave-fraction", type=int, default=12,
                        help="Smoothing bins per octave. Default: 12")
    parser.add_argument("--ylim", type=float, nargs=2, default=None,
                        help="Y-axis limits in dB rel. the reference")
    # Defaults give a ~10.2 x 11 in figure that fills a 6.5 in US-letter
    # text width at about 7 in tall, leaving room for the caption
    parser.add_argument("--panel-width", type=float, default=3.4,
                        help="Panel width in inches. Default: 3.4")
    parser.add_argument("--panel-height", type=float, default=1.7,
                        help="Panel height in inches. Default: 1.7")
    parser.add_argument("--gap", type=float, default=0.01,
                        help="Space between panels as a fraction of the "
                             "panel size. Default: 0.01")
    parser.add_argument("--fontsize", type=float, default=14,
                        help="Base font size in points; sized so text stays "
                             "readable when the figure is shrunk to a 7 in "
                             "page width. Default: 14")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--exclude-sites", default=DEFAULT_EXCLUDED_SITES,
                        help="Comma list of sites to leave out. Pass '' to "
                             f"keep all. Default: {DEFAULT_EXCLUDED_SITES}")
    parser.add_argument("--outfile", default=None,
                        help="Output image. Default: "
                             "FIGURES/permanent_relative.png (plus the "
                             "site list when --sites is used)")
    parser.set_defaults(**defaults)
    args = parser.parse_args()

    def to_key(name):
        """Normalize 'permanent' / suffix input to a station key"""
        name = name.strip()
        return "Permanent" if name.lower() == "permanent" else name.upper()

    args.reference = to_key(args.reference)
    args.omit = {to_key(k) for k in args.omit.split(",") if k.strip()}
    if args.reference in args.omit:
        parser.error("--omit cannot include the --reference station")
    args.sites = set(args.sites.split(",")) if args.sites else None
    args.exclude_sites = {s for s in args.exclude_sites.split(",") if s}
    if args.outfile is None:
        tag = "_" + "_".join(sorted(args.sites)) if args.sites else ""
        args.outfile = f"FIGURES/permanent_relative{tag}.png"
    args.components = [c.strip().upper() for c in args.components.split(",")]
    bad = set(args.components) - set("ZNEH")
    if bad:
        parser.error(f"unknown component(s): {', '.join(sorted(bad))}")
    args.exclude = [c for c in args.exclude.split(",") if c]
    return args


def main(description=__doc__, **defaults):
    """
    Load the spectra and draw the figure.

    :type description: str
    :param description: help text shown by --help
    :type defaults: dict
    :param defaults: overrides for argument defaults (see parse_args)
    """
    args = parse_args(description, **defaults)
    by_site = load_stations(args.spectra, args.sites, args.exclude)
    for site in args.exclude_sites & set(by_site):
        print(f"Leaving out site {site}")
        del by_site[site]
    # Drop omitted stations before anything else so they also don't
    # shorten the common time window
    for site in by_site:
        by_site[site] = [s for s in by_site[site]
                         if station_key(s) not in args.omit]
    if not by_site:
        raise SystemExit("No matching NPZ files found")

    # Rows ordered by the start day of each site's common window, then by
    # name (sites starting on the same day are sorted alphabetically)
    order = sorted(by_site,
                   key=lambda s: (common_window(by_site[s])[0] // 86400, s))
    by_site = {site: by_site[site] for site in order}
    for site, stations in by_site.items():
        names = ", ".join(f"{s.net}.{s.sta}" for s in stations)
        print(f"{site}: {names}")
    plot(by_site, args.components, args.outfile, args)


if __name__ == "__main__":
    main()
