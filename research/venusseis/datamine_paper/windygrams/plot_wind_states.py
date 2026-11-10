#!/usr/bin/env python
"""
Make the wind-state paper figures with PrettyPlot.

All figures are grids where each panel has a spectrogram above the
waveform (response removed, velocity), with the 1-min mean wind speed from
the co-located anemometer on a second axis. Axis limits and spectrogram
colors are shared by every panel of a figure so panels compare directly.

- wind_states.png: K24KU during rising, falling, high and low wind
- wind_states_mitd.png: MITDU during rising, falling, high and low wind
- wind_states_mitd_transients.png: as above, for windows that also contain
  many impulsive transients detected with STA/LTA on the co-located sensors
- wind_states_mitd_transients_1hz.png: as above with all data highpassed
  at 1 Hz, spectrograms starting at 1 Hz, and faint vertical marks on the
  waveforms at the detected transients; wind_states_mitdc_*, _mitdw_* and
  _hvmitd_* show the same for the co-located sensors of MITDU
- wind_compare_mitdu_hvmitd_1hz.png: MITDU (left) next to HV.MITD (right)
  during the rising and high wind windows of those figures
- wind_compare_mitdu_mitdc_hvmitd_1hz.png: as above with MITDC in a middle
  column
- wind_states_stations.png: the same four states, each at a different site
- wind_covers.png, wind_covers_mitd.png: all co-located stations of K24K or
  MITD (Covered, Uncovered, Windscreen, permanent) during rising wind

Data are highpassed at 0.1 Hz, the low end of the band of interest, so the
waveforms show that the noise field below 1 Hz is dominated by the
microseism. Spectrograms use a logarithmic frequency axis starting at
0.5 Hz, which keeps the microseism from washing out the color scale.

The temporary network code of the nodes (XX) is left out of all labels.
Run this script from its symlink in the windygrams working directory.

Example::

    python plot_wind_states.py
    python plot_wind_states.py --component N --fmin 5
"""
import argparse
import csv
import math
import os
from unittest import mock

import matplotlib
matplotlib.use("Agg")  # noqa: E402
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from obspy import UTCDateTime, read

from pysep.pplot import PrettyPlot
from smooth_wind import smooth

PATH_DATA = "DATA"
PATH_WIND = "WIND"
PATH_XML = "XML"
PATH_OUT = "FIGURES/paper"
# Onsets of STA/LTA transients: site, jday, start, offset_s from window start
PATH_EVENTS = "CHECKS/transients/event_times.csv"
TMARK_COLOR = "C1"
TMARK_ALPHA = 0.4
YEAR = 2025
# Node channel code and the permanent station with its channel code
SITES = {
    "DESD": ("HH", "HV.DESD", "EH"),
    "K24K": ("EH", "AK.K24K", "BH"),
    "KAED": ("HH", "HV.KAED", "EH"),
    "MITD": ("HH", "HV.MITD", "HH"),
    "PAX": ("EH", "AK.PAX", "BH"),
    "RND": ("EH", "AK.RND", "BH"),
}
FREQ_TICKS = [0.5, 1, 2, 5, 10, 20]
WIND_YLIM = [0, 14]
WIND_TICKS = [0, 4, 8, 12]

# Panels: (wind state, site, julian day, start, end in UTC)
STATES = [("Rising wind", "K24K", 261, "14:10", "18:10"),
          ("Falling wind", "K24K", 268, "04:50", "08:50"),
          ("High wind", "K24K", 267, "13:10", "17:10"),
          ("Low wind", "K24K", 263, "14:50", "18:50")]
STATES_MITD = [("Rising wind", "MITD", 117, "17:40", "21:40"),
               ("Falling wind", "MITD", 100, "10:30", "14:30"),
               ("High wind", "MITD", 97, "02:40", "06:40"),
               ("Low wind", "MITD", 91, "06:00", "10:00")]
STATES_MITD_TRANSIENTS = [("Rising wind", "MITD", 100, "13:50", "17:50"),
                          ("Falling wind", "MITD", 105, "05:50", "09:50"),
                          ("High wind", "MITD", 95, "07:20", "11:20"),
                          ("Low wind", "MITD", 117, "06:20", "10:20")]
# Output name of the 1 Hz transient figure of each co-located MITD sensor
MARKED_FIGURES = {"U": "wind_states_mitd_transients_1hz.png",
                  "C": "wind_states_mitdc_transients_1hz.png",
                  "W": "wind_states_mitdw_transients_1hz.png",
                  "P": "wind_states_hvmitd_transients_1hz.png"}
STATES_STATIONS = [("Rising wind", "MITD", 117, "17:40", "21:40"),
                   ("Falling wind", "KAED", 104, "00:50", "04:50"),
                   ("High wind", "K24K", 267, "13:10", "17:10"),
                   ("Low wind", "PAX", 247, "06:40", "10:40")]
STATES_STA = "U"
# Rising-wind window of each wind-cover figure: (site, day, start, end)
COVERS_K24K = ("K24K", 261, "14:10", "18:10")
COVERS_MITD = ("MITD", 96, "16:00", "20:00")
COVERS = [("Covered", "C"), ("Uncovered", "U"), ("Windscreen", "W"),
          ("Permanent", "P")]


class PanelPlot(PrettyPlot):
    """
    PrettyPlot that draws into existing axes of a larger figure, leaving
    the figure layout and saving to the caller
    """
    def __init__(self, fids, fig, ax_spectra, ax, wf_fmin=None, **kwargs):
        """
        :type fids: list of str
        :param fids: waveform files, passed to PrettyPlot
        :type fig: matplotlib.figure.Figure
        :param fig: figure holding the axes
        :type ax_spectra: matplotlib.axes.Axes
        :param ax_spectra: axis for the spectrogram
        :type ax: matplotlib.axes.Axes
        :param ax: axis for the waveform
        :type wf_fmin: float or None
        :param wf_fmin: extra highpass corner in Hz for the plotted waveform
            only; the spectrogram keeps the band set by PrettyPlot's `fmin`
        :param kwargs: keyword arguments of PrettyPlot
        """
        super().__init__(fids, **kwargs)
        self._panel = (fig, ax_spectra, ax)
        self.wf_fmin = wf_fmin

    def setup_plot(self):
        """Use the given axes instead of creating a figure"""
        self.f, self.ax_spectra, self.ax = self._panel
        self.axs = [self.ax] * len(self.st)
        self.twax = None

    def plot_waveforms(self):
        """Apply the waveform-only highpass, then plot the waveforms"""
        if self.wf_fmin:
            st = self.st.copy().filter("highpass", freq=self.wf_fmin,
                                       corners=self.corners,
                                       zerophase=self.zerophase)
            self._data = [tr.data for tr in st]
        super().plot_waveforms()

    def plot_tmarks(self):
        """Draw faint time marks on the waveform axis"""
        for tmark in self.tmarks:
            time = mdates.date2num(UTCDateTime(tmark).datetime)
            self.ax.axvline(time, c=TMARK_COLOR, lw=1, alpha=TMARK_ALPHA,
                            zorder=100)

    def finalize(self):
        """Leave the figure open; the caller lays it out and saves it"""


def trace_label(code):
    """
    Return a trace ID without the temporary XX network code

    :type code: str
    :param code: trace ID, e.g. 'XX.K24KU..EHZ'
    :rtype: str
    :return: label, e.g. 'K24KU..EHZ'
    """
    return code[3:] if code.startswith("XX.") else code


def trace_id(site, sta, comp):
    """
    Return the trace ID of one co-located station

    :type site: str
    :param site: site name, key of `SITES`
    :type sta: str
    :param sta: C, U or W for the nodes, P for the permanent station
    :type comp: str
    :param comp: component, Z, N or E
    :rtype: str
    :return: trace ID, e.g. 'XX.K24KU..EHZ'
    """
    node_cha, perm, perm_cha = SITES[site]
    if sta == "P":
        return f"{perm}..{perm_cha}{comp}"
    return f"XX.{site}{sta}..{node_cha}{comp}"


def wind_file(site, jday):
    """
    Return the 1-min mean wind speed file of a day, making it if needed

    :type site: str
    :param site: site name
    :type jday: int
    :param jday: julian day
    :rtype: str
    :return: path to the smoothed wind speed file
    """
    name = f"XX.{site}A..LWS.{YEAR}.{jday:03d}"
    fid = os.path.join(PATH_WIND, name)
    if not os.path.exists(fid):
        smooth(os.path.join(PATH_DATA, name), 60., PATH_WIND)
    return fid


def wind_limit(panels):
    """
    Return the wind speed limit that fits every panel of a figure

    :type panels: list of dict
    :param panels: panels of the figure, see `draw_panel`
    :rtype: int
    :return: peak 1-min mean wind speed in m/s, rounded up
    """
    peak = 0.
    for panel in panels:
        day = f"{YEAR}-{panel['jday']:03d}"
        st = read(wind_file(panel["site"], panel["jday"]),
                  starttime=UTCDateTime(f"{day}T{panel['start']}"),
                  endtime=UTCDateTime(f"{day}T{panel['end']}"))
        peak = max(peak, max(tr.data.max() for tr in st))
    return math.ceil(peak)


def event_times(site, jday, start):
    """
    Return the detected transient onsets of one window

    :type site: str
    :param site: site name
    :type jday: int
    :param jday: julian day
    :type start: str
    :param start: window start as HH:MM in UTC
    :rtype: list of str
    :return: sorted onset times as UTCDateTime strings
    """
    t0 = UTCDateTime(f"{YEAR}-{jday:03d}T{start}")
    with open(PATH_EVENTS, newline="") as f:
        return sorted({str(t0 + float(row["offset_s"]))
                       for row in csv.DictReader(f)
                       if row["site"] == site and int(row["jday"]) == jday
                       and row["start"] == start})


def draw_panel(fig, ax_sp, ax, panel, opts):
    """
    Draw one spectrogram and waveform panel

    :type fig: matplotlib.figure.Figure
    :param fig: figure holding the axes
    :type ax_sp: matplotlib.axes.Axes
    :param ax_sp: spectrogram axis
    :type ax: matplotlib.axes.Axes
    :param ax: waveform axis
    :type panel: dict
    :param panel: 'site', 'code', 'jday', 'start', 'end' and 'title' of
        the panel, optionally 'tmarks' with times to mark
    :type opts: dict
    :param opts: 'fmin', 'wf_fmin', 'ylim', 'sp_vmin', 'sp_vmax', 'sp_ylim'
        and 'sp_log', optionally 'wind_ylim'
    :rtype: tuple
    :return: the wind axis and the colorbar axis
    """
    code, jday = panel["code"], panel["jday"]
    netsta = ".".join(code.split(".")[:2])
    fid_inv = "ISHTAR.xml" if code.startswith("XX.") else f"{netsta}.xml"
    start = UTCDateTime(f"{YEAR}-{jday:03d}T{panel['start']}")
    end = UTCDateTime(f"{YEAR}-{jday:03d}T{panel['end']}")
    before = set(fig.axes)
    pp = PanelPlot(
        [os.path.join(PATH_DATA, f"{code}.{YEAR}.{jday:03d}")],
        fig, ax_sp, ax, wf_fmin=opts["wf_fmin"],
        inv=os.path.join(PATH_XML, fid_inv), merge=True,
        trim_pct=0.1, time="a", rr_output="VEL", fmin=opts["fmin"],
        xlim=[str(start), str(end)], ylim=[-opts["ylim"], opts["ylim"]],
        ylabel="Velocity (m/s)", linewidth=0.3, spectrogram=True,
        sp_vmin=opts["sp_vmin"], sp_vmax=opts["sp_vmax"],
        sp_logscale=opts["sp_log"],
        add_trace=wind_file(panel["site"], jday), tr_ylabel="Wind speed (m/s)",
        tr_color="C3", tr_ylim=opts.get("wind_ylim", WIND_YLIM),
        tr_front=True, legend=False,
        title="", show=False, tmarks=panel.get("tmarks"),
    )
    # Some day files hold other channels too (e.g. AK LHZ next to BHZ)
    pp.st = pp.st.select(channel=code.split(".")[-1])
    nyq = pp.st[0].stats.sampling_rate / 2
    pp.rr_pre_filt = [0.05, 0.1, 0.9 * nyq, nyq]
    # PrettyPlot's tight_layout would rearrange the whole figure each time
    with mock.patch("matplotlib.pyplot.tight_layout"):
        pp.run()

    ax_sp.set_title(panel["title"], loc="left", fontsize=11, pad=3)
    # Grid lines on a spectrogram can be mistaken for spectral lines
    ax_sp.grid(False, which="both")
    if opts["sp_log"]:
        ax_sp.set_yticks([f for f in FREQ_TICKS
                          if opts["sp_ylim"][0] <= f <= opts["sp_ylim"][1]])
        ax_sp.yaxis.set_major_formatter(FuncFormatter(
            lambda v, _: f"{v:g}"))
    ax_sp.set_ylim(opts["sp_ylim"])
    for axis in (ax_sp, ax):
        axis.xaxis.set_major_locator(mdates.HourLocator())
        axis.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    new = [a for a in fig.axes if a not in before and a.axison]
    cax = [a for a in new if a is not pp.twax][0]
    return pp.twax, cax


def style_panel(ax_sp, ax, twax, cax, left, right, bottom, ylim, yticks=None,
                wind_ticks=None, xlabel=None):
    """
    Set shared-axis labels and font sizes of one panel for print

    :type ax_sp: matplotlib.axes.Axes
    :param ax_sp: spectrogram axis
    :type ax: matplotlib.axes.Axes
    :param ax: waveform axis
    :type twax: matplotlib.axes.Axes
    :param twax: wind axis
    :type cax: matplotlib.axes.Axes
    :param cax: colorbar axis
    :type left: bool
    :param left: panel is in the left column
    :type right: bool
    :param right: panel is in the right column
    :type bottom: bool
    :param bottom: panel is in the bottom row
    :type ylim: float
    :param ylim: symmetric waveform y-limit in m/s, sets the axis scale
    :type yticks: list of float or None
    :param yticks: waveform y-ticks in m/s, automatic if None
    :type wind_ticks: list of int or None
    :param wind_ticks: wind speed ticks in m/s, `WIND_TICKS` if None
    :type xlabel: str or None
    :param xlabel: time axis label for this panel in any row, if None the
        bottom row is labeled 'Time (UTC)'
    """
    for axis in (ax_sp, ax, twax, cax):
        axis.tick_params(labelsize=10, length=3.5, width=0.8)
        axis.tick_params(which="minor", length=1.5, width=0.5)
        for spine in axis.spines.values():
            spine.set_linewidth(0.8)
    ax_sp.tick_params(axis="x", labelbottom=False)
    ax_sp.set_xlabel("")
    ax_sp.set_ylabel("Freq. (Hz)" if left else "", fontsize=11)
    # The scale goes in the label: an offset text would hit the spectrogram
    exp = math.floor(math.log10(ylim))
    if yticks:
        ax.set_yticks(yticks)
    ax.yaxis.set_major_formatter(FuncFormatter(
        lambda v, _: f"{v * 10. ** -exp:g}".replace("-", "\u2212")))
    ax.set_ylabel(f"Velocity (10$^{{{exp}}}$ m/s)" if left else "",
                  fontsize=11)
    ax_sp.tick_params(axis="y", labelleft=left)
    ax.tick_params(axis="y", labelleft=left)
    ax.set_xlabel(xlabel or ("Time (UTC)" if bottom else ""), fontsize=11)

    twax.set_ylabel("Wind speed (m/s)" if right else "", fontsize=11,
                    color="C3", rotation=270,
                    labelpad=11 if right else 0)
    twax.set_yticks(wind_ticks or WIND_TICKS)
    twax.tick_params(axis="y", labelright=right, colors="C3")
    cax.tick_params(axis="y", labelright=right)
    cax.set_ylabel(r"PSD (dB re (m/s)$^2$/Hz)" if right else "", fontsize=10,
                   rotation=270, labelpad=12 if right else 0)


def make_figure(panels, fid_out, opts):
    """
    Draw a grid of spectrogram and waveform panels and save it

    :type panels: list of dict
    :param panels: the panels in reading order, see `draw_panel`
    :type fid_out: str
    :param fid_out: output figure path
    :type opts: dict
    :param opts: plotting options, see `draw_panel`; 'wind_auto' scales the
        wind axis to the highest wind speed of the figure, 'ncols' sets
        the number of columns (default 2), 'hspace' the space between rows
        (default 0.2) and 'date_xlabel' adds YYYY-DDD to every time label
    :rtype: str
    :return: path of the saved figure
    """
    opts = dict(opts)
    if opts.get("wind_auto"):
        top = wind_limit(panels)
        opts["wind_ylim"] = [0, top]
        opts["wind_ticks"] = list(range(0, top + 1, 2))
    # 3:2 aspect ratio; fonts are sized for reduction to the page width
    fig = plt.figure(figsize=(10.5, 7.), dpi=300)
    ncols = opts.get("ncols", 2)
    nrows = math.ceil(len(panels) / ncols)
    outer = fig.add_gridspec(nrows, ncols, left=0.065, right=0.905, top=0.955,
                             bottom=0.08, hspace=opts.get("hspace", 0.2),
                             wspace=0.07)
    for i, panel in enumerate(panels):
        inner = outer[i].subgridspec(2, 1, hspace=0.06)
        ax_sp = fig.add_subplot(inner[0])
        ax = fig.add_subplot(inner[1], sharex=ax_sp)
        twax, cax = draw_panel(fig, ax_sp, ax, panel, opts)
        style_panel(ax_sp, ax, twax, cax, left=i % ncols == 0,
                    right=i % ncols == ncols - 1,
                    bottom=i >= (nrows - 1) * ncols, ylim=opts["ylim"],
                    yticks=opts.get("yticks"),
                    wind_ticks=opts.get("wind_ticks"),
                    xlabel=(f"Time (UTC), {YEAR}-{panel['jday']:03d}"
                            if opts.get("date_xlabel") else None))
    os.makedirs(os.path.dirname(fid_out), exist_ok=True)
    fig.savefig(fid_out, dpi=300)
    plt.close(fig)
    return fid_out


def date_label(jday):
    """
    Return the calendar date of a julian day

    :type jday: int
    :param jday: julian day
    :rtype: str
    :return: date as YYYY-MM-DD
    """
    return str(UTCDateTime(f"{YEAR}-{jday:03d}"))[:10]


def states_panels(states, sta, comp, marks=False):
    """
    Return the panels of the wind-state figure

    :type states: list of tuple
    :param states: (wind state, site, julian day, start, end) of each panel
    :type sta: str
    :param sta: station, see `trace_id`
    :type comp: str
    :param comp: component
    :type marks: bool
    :param marks: mark the detected transients listed in `PATH_EVENTS`
    :rtype: list of dict
    :return: panels for `make_figure`
    """
    panels = []
    for i, (state, site, jday, t0, t1) in enumerate(states):
        code = trace_id(site, sta, comp)
        panels.append({"site": site, "code": code, "jday": jday, "start": t0,
                       "end": t1,
                       "tmarks": (event_times(site, jday, t0) if marks
                                  else None),
                       "title": f"({chr(97 + i)}) {state}, "
                                f"{trace_label(code)}, {date_label(jday)}"})
    return panels


def compare_panels(states, stations, comp, dates=True):
    """
    Return panels with one wind state per row and one station per column

    :type states: list of tuple
    :param states: (wind state, site, julian day, start, end) of each row
    :type stations: list of str
    :param stations: station of each column, see `trace_id`
    :type comp: str
    :param comp: component
    :type dates: bool
    :param dates: include the date in panel titles
    :rtype: list of dict
    :return: panels for `make_figure`, with detected transients marked
    """
    panels = []
    for state, site, jday, t0, t1 in states:
        for sta in stations:
            code = trace_id(site, sta, comp)
            title = (f"({chr(97 + len(panels))}) {state}, "
                     f"{trace_label(code)}")
            panels.append({"site": site, "code": code, "jday": jday,
                           "start": t0, "end": t1,
                           "tmarks": event_times(site, jday, t0),
                           "title": (f"{title}, {date_label(jday)}" if dates
                                     else title)})
    return panels


def covers_panels(window, comp):
    """
    Return the panels of the wind-cover figure

    :type window: tuple
    :param window: (site, julian day, start, end) of the rising-wind window
    :type comp: str
    :param comp: component
    :rtype: list of dict
    :return: panels for `make_figure`
    """
    site, jday, t0, t1 = window
    panels = []
    for i, (cover, sta) in enumerate(COVERS):
        code = trace_id(site, sta, comp)
        panels.append({"site": site, "code": code, "jday": jday,
                       "start": t0, "end": t1,
                       "title": f"({chr(97 + i)}) {cover}, "
                                f"{trace_label(code)}"})
    return panels


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--component", default="Z", choices=["Z", "N", "E"])
    parser.add_argument("--fmin", type=float, default=0.1,
                        help="highpass corner in Hz, default 0.1")
    parser.add_argument("--wf_fmin", type=float, default=None,
                        help="extra highpass corner in Hz of the plotted "
                             "waveform only, default none")
    parser.add_argument("--sp_fmin", type=float, default=0.5,
                        help="lowest spectrogram frequency in Hz, default "
                             "0.5")
    parser.add_argument("--sp_linear", action="store_true",
                        help="use a linear spectrogram frequency axis")
    parser.add_argument("--ylim_states", type=float, default=4E-5,
                        help="waveform y-limit in m/s, wind-state figure")
    parser.add_argument("--ylim_mitd", type=float, default=2E-5,
                        help="waveform y-limit in m/s, MITDU wind-state "
                             "figure")
    parser.add_argument("--ylim_marked", type=float, default=5E-6,
                        help="waveform y-limit in m/s, 1 Hz MITDU figure "
                             "with transient marks")
    parser.add_argument("--ylim_covers", type=float, default=2E-5,
                        help="waveform y-limit in m/s, wind-cover figure")
    parser.add_argument("--sp_vmin", type=float, default=-190.)
    parser.add_argument("--sp_vmax", type=float, default=-110.)
    parser.add_argument("--path_out", default=PATH_OUT)
    return parser.parse_args()


def main():
    """Make all of the wind-state and wind-cover figures"""
    args = parse_args()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    matplotlib.rcParams["agg.path.chunksize"] = 10000
    opts = {"fmin": args.fmin, "wf_fmin": args.wf_fmin,
            "ylim": args.ylim_states,
            "sp_vmin": args.sp_vmin, "sp_vmax": args.sp_vmax,
            "sp_ylim": [args.sp_fmin, 45.], "sp_log": not args.sp_linear}
    comp = args.component
    for fid, states, ylim in (
            ("wind_states.png", STATES, args.ylim_states),
            ("wind_states_mitd.png", STATES_MITD, args.ylim_mitd),
            ("wind_states_mitd_transients.png", STATES_MITD_TRANSIENTS,
             args.ylim_mitd),
            ("wind_states_stations.png", STATES_STATIONS, args.ylim_states)):
        opts["ylim"] = ylim
        print(make_figure(states_panels(states, STATES_STA, comp),
                          os.path.join(args.path_out, fid), opts))
    marked = {**opts, "fmin": 1., "sp_ylim": [1., 45.],
              "ylim": args.ylim_marked, "yticks": [-4E-6, 0., 4E-6],
              "wind_auto": True}
    for sta, fid in MARKED_FIGURES.items():
        print(make_figure(
            states_panels(STATES_MITD_TRANSIENTS, sta, comp, marks=True),
            os.path.join(args.path_out, fid), marked))
    # Rising and high wind only, uncovered node next to the permanent station
    print(make_figure(
        compare_panels([STATES_MITD_TRANSIENTS[0], STATES_MITD_TRANSIENTS[2]],
                       ["U", "P"], comp),
        os.path.join(args.path_out, "wind_compare_mitdu_hvmitd_1hz.png"),
        marked))
    print(make_figure(
        compare_panels([STATES_MITD_TRANSIENTS[0], STATES_MITD_TRANSIENTS[2]],
                       ["U", "C", "P"], comp, dates=False),
        os.path.join(args.path_out,
                     "wind_compare_mitdu_mitdc_hvmitd_1hz.png"),
        {**marked, "ncols": 3, "date_xlabel": True, "hspace": 0.28}))
    opts["ylim"] = args.ylim_covers
    print(make_figure(covers_panels(COVERS_MITD, comp),
                      os.path.join(args.path_out, "wind_covers_mitd.png"),
                      opts))
    # AK.K24K records up to 25 Hz, so all panels of its figure stop there
    opts["sp_ylim"] = [args.sp_fmin, 22.5]
    print(make_figure(covers_panels(COVERS_K24K, comp),
                      os.path.join(args.path_out, "wind_covers.png"), opts))


if __name__ == "__main__":
    main()
