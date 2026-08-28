"""
Calculate and scatter plot P-to-S amplitude ratio for a set of synthetics.
Timing windows based on TauP. Allows inputting multiple sets of synthetics
to compare different events etc.
"""
import argparse
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from matplotlib import patheffects as pe
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from obspy.imaging.beachball import beach
from obspy.taup import TauPyModel
from obspy.geodetics import kilometers2degrees, gps2dist_azimuth
from pysep import read_sem, logger

logger.setLevel("CRITICAL")


ARRIVAL_RG=3.1
FID_COAST = "/home/bhchow/work/data/cartography/coastline_128_130_40_43.csv"

def get_p2s(tr, p_window, s_window, choice="window"):
    """
    Calculate P-to-S amplitude ratio for one synthetic
    
    :param: ObsPy Trace holding data
    :param p_window: (start, stop) in units of s to search for P-arrival
    :param s_window: (start, stop) in units of s to search for S-arrival
    :param fmin: min freq. for filter in [Hz]
    :param fmax: max freq. for filter in [Hz]
    :params choice: 
        - 'window' for looking at start and stop of phase list (conservative)
        - 'bounds' for searching before latest P and after earliest S (loose)
    :return: absolute amplitude ratio
    """
    sr = tr.stats.sampling_rate
    p_start, p_stop = [_ * sr for _ in p_window]  # sec -> Hz
    s_start, s_stop = [_ * sr for _ in s_window]
    
    p_start = int(p_start)
    p_stop = int(p_stop)
    if p_stop == p_start:
        p_stop += 1
    s_start = int(s_start)
    s_stop = int(s_stop)
    if s_stop == s_start:
        s_stop += 1

    # Only look inside the given P or S window
    try:
        if choice == "window":
            p_idx = np.argmax(tr[p_start: p_stop]) + p_start
            s_idx = np.argmax(tr[s_start: s_stop]) + s_start
        # Use the extremum bounds to search for P and S
        elif choice == "bounds":
            # Slices start at sample 0, so argmax already returns an absolute
            # index -- do NOT add p_start/s_start here (that was double-offsetting)
            p_idx = np.argmax(tr[:p_stop])
            s_idx = np.argmax(tr[s_start:]) + s_start
        # Use min S time as a dividing line
        elif choice == "before_after_s":
            p_idx = np.argmax(tr[:s_start])
            s_idx = np.argmax(tr[s_start:]) + s_start
    except ValueError:
        breakpoint()

    # P-to-S amplitude ratio, index of max P, index of max S
    return float(np.abs(tr[p_idx] / tr[s_idx])), p_idx, s_idx


def get_taup_arrivals(source_depth_in_km, distance_in_km, p_phase_list=None,
                      s_phase_list=None, 
                      model="/Users/prof/Repos/spectral/research/seismology/"
                            "taup_models/ak135f_upper_crust.npz"):
    """
    Get arrival time windows from TauP for a given `TAUP_MODEL`
    """
    # By default query all the crustal directarrivals
    if not p_phase_list:
        p_phase_list = ["p", "P"]
    if not s_phase_list:
        s_phase_list = ["s", "S"]

    dist_deg = kilometers2degrees(distance_in_km)

    model = TauPyModel(model=model)
    p_arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km,
                                        distance_in_degree=dist_deg,
                                        phase_list=p_phase_list)

    s_arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km,
                                        distance_in_degree=dist_deg,
                                        phase_list=s_phase_list)

    if not p_arrivals or not s_arrivals:
        return None, None
    
    # Aggregate all P and S-wave arrivals
    p_arrival_times = [_.time for _ in p_arrivals]
    s_arrival_times = [_.time for _ in s_arrivals]

    # Window from the start of the P-wave to the start of the S-wave
    p_window = [min(p_arrival_times), min(s_arrival_times)]

    # Window from start of S-wave to end of S-wave
    # Set a backstop velocity so that we don't pick into Rg @ 3.2km/s
    s_window = [min(s_arrival_times), max(s_arrival_times)]
    arrival_rg = distance_in_km / ARRIVAL_RG  # s
    if s_window[1] > arrival_rg:
        s_window[1] = arrival_rg
    if s_window[1] < s_window[0]:
        s_window[1] = s_window[0] 

    return p_window, s_window, p_arrivals, s_arrivals


def plot_tr(tr, p_idx=None, s_idx=None, p_window=None, s_window=None,
            p_arrivals=None, s_arrivals=None,
            title_prepend="", save="./", show=False):
    """Make individual waveform plot"""
    # FOR WAVEFORM PLOTS, NOT THE HEATMAP
    f, ax = plt.subplots(1, figsize=(8,6), dpi=100)

    plt.plot(tr.times(), tr.data, c="k", zorder=9, lw=1)

    if p_idx:
        plt.scatter(tr.times()[p_idx], tr.data[p_idx], zorder=10, c="C0", 
                    label=f"P_max={tr.data[p_idx]:.2E}; "
                          f"t={tr.times()[p_idx]:.2f}")
    if s_idx:
        plt.scatter(tr.times()[s_idx], tr.data[s_idx], zorder=10, c="C1", 
                    label=f"S_max={tr.data[s_idx]:.2E}; "
                          f"t={tr.times()[s_idx]:.2f}")

    if p_arrivals:
        for arrival in p_arrivals:
            plt.axvline(arrival.time, c="C0", alpha=0.5,)
            plt.text(arrival.time, tr.data.max() / 2, s=arrival.name)

    if s_arrivals:
        for arrival in s_arrivals:
            plt.axvline(arrival.time, c="C1", alpha=0.5,)
            plt.text(arrival.time, tr.data.max() / 2, s=arrival.name)
    

    if p_window:
        p_patch = Rectangle(xy=(p_window[0], tr.data.min()),
                            width=p_window[1]-p_window[0], 
                            height=tr.data.max() - tr.data.min(),
                            color="C0", alpha=0.25)
        ax.add_patch(p_patch)
    if s_window:
        s_patch = Rectangle(xy=(s_window[0], tr.data.min()),
                            width=s_window[1]-s_window[0], 
                            height=tr.data.max() - tr.data.min(),
                            color="C1", alpha=0.25)
        ax.add_patch(s_patch)

    if s_window:
        xmax = min(s_window[1] + 5, tr.times()[-1])
    else:
        xmax = tr.times()[-1]

    if p_window:
        xmin = max(p_window[0] - 5, 0)
    else:
        xmin = 0
    plt.xlim(xmin, xmax)

    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel("Velocity [m/s]")
    dist = tr.stats.sac["dist"]
    baz = tr.stats.sac["baz"]
    plt.title(f"{title_prepend}; "
              f"{tr.get_id()}; "
              f"Dist={dist:.2f}km; "
              f"BAz={baz:.2f}deg\n"
              f"P/S={tr.data[p_idx] / tr.data[s_idx]:.2f}"
              )

    if show:
        plt.show()
    if save:
        plt.savefig(os.path.join(save, f"{tr.get_id()}.png"))
    plt.close(f)


def add_waveform_inset(fig, ax, lon, lat, tr, p_idx=None, s_idx=None,
                       width=0.11, height=0.075, pad_s=1.0,
                       annotate=True):
    """
    Draw a small, frameless waveform (with P/S picks marked) whose visible
    window is pinned exactly at a station's map location -- just the trace,
    pick markers, and two small annotations floating directly on top of the
    heatmap, no box/background/title/leader line.

    The inset is cropped to the P-to-S window (+/- `pad_s`), and both axes
    are scaled to that visible slice of data rather than the full trace, so
    a long, mostly-flat trace doesn't shrink the interesting part down to a
    sliver.

    :param fig: parent matplotlib Figure
    :param ax: parent map Axes (in lon/lat data coordinates)
    :param lon/lat: station location the waveform's visible start is pinned to
    :param tr: ObsPy Trace to plot
    :param p_idx/s_idx: sample indices of the picked P/S peak amplitudes
    :param width/height: inset size in figure-fraction units
    :param pad_s: seconds of padding before the P pick / after the S pick
    """
    times = tr.times()
    data = tr.data

    xlo = max(times[0], times[p_idx] - pad_s) if p_idx is not None else times[0]
    # xhi = min(times[-1], times[s_idx] + pad_s) if s_idx is not None else times[-1]

    # Stop after a given group velocity
    xhi =  tr.stats.sac.dist / 3.0

    # Scale to only the data actually shown in [xlo, xhi], not the full trace
    mask = (times >= xlo) & (times <= xhi)
    window_data = data[mask]
    dmin, dmax = window_data.min(), window_data.max()
    pad = 0.1 * ((dmax - dmin) or 1)
    ylo, yhi = dmin - pad, dmax + pad

    # Anchor point: the trace's (interpolated) value at the left edge of the
    # visible window -- this is what gets pinned to the station marker
    y_at_xlo = np.interp(xlo, times, data)
    y_frac = (y_at_xlo - ylo) / (yhi - ylo)

    # Station location -> figure-fraction coordinates, so the inset can be
    # placed as a sibling Axes rather than embedded inside the map Axes
    fx, fy = fig.transFigure.inverted().transform(
        ax.transData.transform((lon, lat)))

    # Anchor the inset so (xlo, y_at_xlo) lands exactly on the marker: left
    # edge at the marker's x, and shifted vertically so y_frac lines up
    x0 = min(max(fx, 0.01), 0.99 - width)
    y0 = min(max(fy - y_frac * height, 0.01), 0.99 - height)

    inset_ax = fig.add_axes([x0, y0, width, height])
    inset_ax.plot(times, data, c="lime", lw=1.5, zorder=1)

    if p_idx is not None:
        inset_ax.scatter(times[p_idx], data[p_idx], c="y", s=35, zorder=2,
                         ec="k", marker="X")
    if s_idx is not None:
        inset_ax.scatter(times[s_idx], data[s_idx], c="c", s=35, zorder=2,
                         marker="X", ec="k")

    inset_ax.set_xlim(xlo, xhi)
    inset_ax.set_ylim(ylo, yhi)

    # No frame: no ticks, no spines, no background patch -- just the trace
    # and P/S markers floating directly on top of the heatmap
    inset_ax.axis("off")
    inset_ax.patch.set_alpha(0)

    # Annotate the two quantities the inset's axes would otherwise convey:
    # the time span shown, and the peak amplitude shown
    if annotate:
        bbox = dict(boxstyle="round,pad=0.2", facecolor="w", edgecolor="None",
                    linewidth=0.6, alpha=1)
        inset_ax.annotate(f"{xhi - xlo:.1f}s\n{tr.stats.sac.dist:.2f}km ",
                          xy=(-.02, 0.5),
                          xycoords="axes fraction", ha="right", va="bottom",
                          fontsize=14, color="lime",
                          path_effects=[pe.withStroke(linewidth=2, foreground="k")]
                          )# bbox=bbox, zorder=3)

    return inset_ax


def add_ruler(ax, point1, point2, color="lime", lw=1.5, fontsize=18):
    """
    Draw a line between two (lon, lat) points on a map Axes and annotate the
    great-circle distance between them in km.

    :param ax: map Axes (in lon/lat data coordinates)
    :param point1/point2: (lon, lat) tuples marking the ruler's endpoints
    :return: distance between the two points, in km
    """
    lon1, lat1 = point1
    lon2, lat2 = point2
    dist_m, _, _ = gps2dist_azimuth(lat1, lon1, lat2, lon2)
    dist_km = dist_m / 1000.

    ax.plot([lon1, lon2], [lat1, lat2], c=color, lw=lw, zorder=15,
           marker="o", markersize=4, markerfacecolor=color,
           markeredgecolor="w")

    ax.annotate(f"{dist_km:.1f} km", xy=((lon1 + lon2) / 2, (lat1 + lat2) / 2),
               xycoords="data", ha="center", va="bottom", fontsize=fontsize,
               color=color, zorder=16)

    return dist_km


class P2SRatio:
    """Gather P2S ratio in class attributes"""
    def __init__(self, path, path_src="CMTSOLUTION", path_sta="STATIONS", 
                 choice="taup", path_save="data.npz", path_fig="figures", 
                 fmin=1, fmax=6, components="ZNE", overwrite=False):
        """
        Setup attributes for P2S

        .. variables::
        path (str): path to waveform data
            path_src (str): path to CMTSOLUTION for finding az/dist
            path_sta (str): path to STATIONS for getting lat/lon
            path_save (str): path and filename for .npz file for saving state
            path_fig: (str): path for saving output figures
            fmin (float): minimum freq corner
            fmax (float): maximum freq corner
            components (str): which components to include, 'Z' for vertical only
                'ZNE' to average all 3.
            overwrite (bool): if True, overwrite data in `path_data`
        """
        self.path = path
        self.tag = str(self.path.parent) + "_" + self.path.name 
        self.choice = choice

        self.path_save = Path(path_save)
        self.path_src = Path(path_src)
        self.path_sta = Path(path_sta)
        self.path_fig = Path(path_fig)

        assert(os.path.exists(self.path_sta)), f"{self.path_sta} does not exist"
        assert(os.path.exists(self.path_src)), f"{self.path_src} does not exist"
        os.makedirs(self.path_fig, exist_ok=True)

        self.fmin = fmin
        self.fmax = fmax
        self.components = components
        self.overwrite = overwrite

        # Only use the first component to get IDs, even if we have multi-comp.
        self.fids = sorted(self.path.glob(f"*.?X[{components[0]}].sem?"))
        assert(self.fids), f"no files found for {self.path}"

        # Get moment tensor components
        self.cmt = self._get_mt()

        # Internal parameter lists
        self.p2sratios = []
        self.lats = []
        self.lons = []
        self.distances = []
        self.ids = []
        self.srclat = None
        self.srclon = None

    def write(self):        
        """Write to .npz file"""
        dict_out = {"p2sratios": np.array(self.p2sratios),
                    "lats": np.array(self.lats),
                    "lons": np.array(self.lons),
                    "distances": np.array(self.distances),
                    "ids": np.array(self.ids),
                    }
        np.savez(self.path_save, **dict_out)
    
    def read(self):
        """Read from .npz file"""
        print(f"loading data from file `{self.path_save}`")
        dict_in = np.load(self.path_save)
        self.p2sratios = dict_in["p2sratios"]
        self.lats = dict_in["lats"]
        self.lons = dict_in["lons"]
        self.distances = dict_in["distances"]   
        self.ids = dict_in["ids"] 

    def _read_sem(self, fid):
        """Wrapper to read sem without calling all args"""
        tr = read_sem(fid, source=self.path_src, stations=self.path_sta)[0]
        # Trim off the T0 values to avoid any confusion
        starttime = tr.stats.starttime - tr.stats.sac.b
        endtime = tr.stats.endtime 
        tr.trim(starttime, endtime)
        return tr

    def _get_mt(self):
        """Get moment tensor components from CMTSOLUTION for plotting"""
        with open(self.path_src) as f:
            lines = f.readlines()
        dict_out = {}
        for line in lines[1:]:
            name, val = line.strip().split(":")
            try:
                dict_out[name] = float(val)
            except ValueError:
                dict_out[name] = val.strip()
        return dict_out
    
    def calculate_ratio(self, i=0, j=-1, parallel=True, ntasks=None):
        """
        Read or fetch P2S amplitude ratios and meta data from the data files

        Args
            i (int): optional, starting index for fids
            j (int): optional, ending index for fids
            parallel (bool): if True, use concurrent futures to speed up process
        """
        if self.path_save.exists() and (not self.overwrite):
            print("reading from database")
            self.read()
        else:
            if not parallel:
                print(f"calculating p2s ratios in serial for {len(self.fids)}")
                error = 0
                for fid in self.fids[i:j]:
                    tr, p2sratio = self.calculate_ratio_single(fid)
                    self.p2sratios.append(p2sratio)
                    self.ids.append(tr.stats.station)
                    self.lats.append(tr.stats.sac["stla"])
                    self.lons.append(tr.stats.sac["stlo"])
                if error:
                    print(f"{self.tag}: {error}/{len(self.fids[i:j])} errors")
            # Use concurrent.futures to parallelize processing
            else:
                print(f"calculating p2s ratios in parallel for {len(self.fids)}")
                with ProcessPoolExecutor(
                    max_workers=ntasks or os.cpu_count()) as executor:
                    futures = [
                        executor.submit(self.calculate_ratio_single, fid) for 
                            fid in self.fids[i:j]
                            ]
                # Collect results
                for future in as_completed(futures):
                    tr, p2sratio = future.result()
                    self.p2sratios.append(p2sratio)
                    self.ids.append(tr.stats.station)
                    self.lats.append(tr.stats.sac["stla"])
                    self.lons.append(tr.stats.sac["stlo"])

            self.write()

        # Last minute, grab source information
        tr = self._read_sem(self.fids[0])
        self.srclat = tr.stats.sac["evla"]
        self.srclon = tr.stats.sac["evlo"]

    def calculate_ratio_single(self, fid):
        """
        Read files and collect information including P2S for a single station.
        `Components` control how many waveforms are opened and included per 
        station

        Args
            fid (str): file identifier pointing to one of the synthetic 
                waveforms for use in analysis
        """
        _corners = 4
        tr = self._read_sem(fid)
        tr.filter("bandpass", freqmin=self.fmin, freqmax=self.fmax,
                  corners=_corners)

        p_window, s_window, p_arrivals, s_arrivals = \
                get_taup_arrivals(tr.stats.sac["evdp"], tr.stats.sac["dist"])

        title_prepend = f"{self.tag} {self.fmin}-{self.fmax}Hz "

        # If no windows can be picked, plot and return
        if not p_window or not s_window:
            if self.path_fig:
                plot_tr(tr, title_prepend=title_prepend, save=self.path_fig)
            print(f"no valid P/S window for {fid}")
            return tr, np.nan
        
        p2sratio, p_idx, s_idx = get_p2s(tr=tr, p_window=p_window, 
                                         s_window=s_window)  
        if self.path_fig:
            plot_tr(tr, p_idx, s_idx, p_window, s_window, 
                    p_arrivals, s_arrivals,
                    title_prepend=title_prepend,
                    save=self.path_fig, show=False)

        # If more than 1 component we need to get the other 2 arrays
        if len(self.components) > 1:
            comp = self.components[0]
            for comp_ in self.components[1:]:
                fid_ = str(fid).replace(f"X{comp}.sem", f"X{comp_}.sem")
                tr_ = self._read_sem(fid_)
                tr_.filter("bandpass", freqmin=self.fmin, freqmax=self.fmax,
                           corners=_corners)

                p2sratio_, p_idx_, s_idx_ = get_p2s(tr=tr_, p_window=p_window, 
                                                    s_window=s_window)  
                p2sratio += p2sratio_
                if self.path_fig:
                    plot_tr(tr_, p_idx_, s_idx_, p_window, s_window,
                            p_arrivals, s_arrivals,
                            title_prepend=title_prepend,
                            save=self.path_fig)
            
            p2sratio /= len(self.components)  # take average

        return tr, p2sratio

    def get_station_waveform(self, station):
        """
        Re-read and process the primary-component waveform for a single
        named station, for overlaying its trace + P/S picks on the heatmap.

        :param station: station code (e.g., 'BFZ') matching the STATIONS file
        :return: (tr, p_idx, s_idx, p_window, s_window) or None if the
            station/windows can't be found
        """
        matches = [fid for fid in self.fids
                   if fid.name.split(".")[1] == station]
        if not matches:
            print(f"no waveform found for station '{station}', skipping")
            return None

        tr = self._read_sem(matches[0])
        tr.filter("bandpass", freqmin=self.fmin, freqmax=self.fmax,
                  corners=4)

        p_window, s_window, *_ = get_taup_arrivals(tr.stats.sac["evdp"],
                                                   tr.stats.sac["dist"])

        if not p_window or not s_window:
            print(f"could not determine windows for station '{station}', "
                  f"skipping")
            return None

        _, p_idx, s_idx = get_p2s(tr=tr, p_window=p_window, s_window=s_window)

        return tr, p_idx, s_idx, p_window, s_window


def plot_scatterplot(paths, fmin=2, fmax=4, components="Z", j=-1):
    """
    Run through P2S for multiple events and create a scatterplot of distance
    versus P2S amplitude ratio. Similar to Walter et al 2018 F2a but without
    magnitude on the x axis
    """
    ratios = []
    for path in paths:
        p2s = P2SRatio(path, path_src=f"CMTSOLUTION_{Path(path).name}",
                       fmin=fmin, fmax=fmax, components=components)
        print(p2s.tag)
        p2s.calculate_ratio(j=j)
        ratios.append(p2s)

    # Generate Scatterplot
    for i, ratio in enumerate(ratios):
        if ratio.tag.startswith("EQ"):
            marker = "o"
            size = 20
        else:
            marker = "*"
            size = 40
        plt.scatter(ratio.distances, ratio.p2sratios, marker=marker, 
                    color=f"C{i}", s=size, label=ratio.tag)

    plt.legend()
    plt.title(f"P2S Ratio [{fmin}-{fmax}]Hz")
    plt.yscale("log")
    plt.xlabel("Distance [km]")
    plt.ylabel("P-to-S Amplitude Ratio")


def plot_heatmap(p2s, threshold=0.8, interpolate=True, method="linear",
                 grid_res=200, smooth_sigma=1.5, save="./figures",
                 cmap="viridis", mt_color="orange", coast_color="k",
                 fid_coastline=FID_COAST, highlight_stations=None,
                 ruler=None, title=None, subplot_label=None, annotate=True,
                 colorbar=True, show=True, log_scale=False):
    """
    For a single event plot a map of amplitude ratios for each station
    interpolated to create a continuous figure rather than a scatterplot.

    :param grid_res: number of interpolation grid points per axis. Fixed
        (not tied to station count) so the heatmap smoothness doesn't
        degrade for sparse/irregular station networks.
    :param smooth_sigma: sigma (in grid cells) of the Gaussian filter applied
        to the interpolated grid to remove small-scale interpolation
        artifacts (e.g. 'cubic' ringing, hard convex-hull edges). Set to 0
        to disable.
    :param highlight_stations: optional list of station codes to annotate on
        the map with a marker + inset plot of their waveform and P/S picks
    :param ruler: optional (point1, point2) pair, or list of such pairs,
        where each point is an (lon, lat) tuple. Draws a line between the
        two points and annotates the great-circle distance between them
    :param log_scale: if True, color the heatmap on a log scale instead of
        linear. P/S ratios often span more than a decade, so log scale can
        resolve structure in the small-ratio (earthquake-like) end that a
        linear scale compresses down near zero.
    """
    print("plotting heatmap")
    f, ax = plt.subplots(1, figsize=(12, 13), dpi=100)

    # Interpolation causes some NaN errors with large values
    # Create a regular grid from lon/lat ranges
    # Set colorscale bounds based on user-defined threshold
    if log_scale:
        # 0 is undefined on a log scale -- use a small positive floor
        # instead, one decade below `threshold`
        vmin = threshold / 10
        vmax = threshold * 2
        levels = np.geomspace(vmin, vmax, int((vmax - vmin) * 10 * 6))
    else:
        vmin = 0
        vmax = vmin + threshold * 2
        levels = np.linspace(vmin, vmax, int((vmax - vmin) * 10 * 6))
    ticks = [vmin, threshold, vmax]

    # Fixed, dense grid independent of station density -- interpolating
    # onto a grid only ~4x the station spacing (the old approach) leaves
    # contourf with little to work with and produces a faceted/jagged
    # look, especially for sparse or irregularly-spaced networks
    pad_lon = 0.05 * (max(p2s.lons) - min(p2s.lons) or 1)
    pad_lat = 0.05 * (max(p2s.lats) - min(p2s.lats) or 1)
    lon_grid = np.linspace(min(p2s.lons) - pad_lon,
                           max(p2s.lons) + pad_lon, grid_res)
    lat_grid = np.linspace(min(p2s.lats) - pad_lat,
                           max(p2s.lats) + pad_lat, grid_res)
    x, y = np.meshgrid(lon_grid, lat_grid)

    # Interpolate within the station convex hull...
    points = (p2s.lons, p2s.lats)
    z_interior = griddata(points, p2s.p2sratios, (x, y), method=method)
    # ...and fill outside it with nearest-neighbor extrapolation instead
    # of a flat fill_value, which used to draw a hard, jagged-looking
    # edge that traced the (irregular) convex hull of the station network
    z_exterior = griddata(points, p2s.p2sratios, (x, y), method="nearest")
    z = np.where(np.isnan(z_interior), z_exterior, z_interior)

    # Light smoothing to knock down remaining small-scale interpolation
    # artifacts without washing out real spatial structure
    if smooth_sigma:
        z = gaussian_filter(z, sigma=smooth_sigma)

    # Colorbar that respects the vmin vmax values. On a log scale, data
    # below `vmin` is much more likely (earthquake-like ratios can sit well
    # under threshold/10), so check both ends rather than just the top
    extend_max = vmax < z.max()
    extend_min = log_scale and (vmin > z.min())
    if extend_max and extend_min:
        extend = "both"
    elif extend_max:
        extend = "max"
    elif extend_min:
        extend = "min"
    else:
        extend = "neither"

    # PLOT
    if log_scale:
        # LogNorm requires strictly positive data -- clip instead of letting
        # any zero/negative values (e.g. from interpolation/smoothing) error
        z_plot = np.clip(z, vmin * 1e-2, None)
        cf = plt.contourf(x, y, z_plot, levels=levels, cmap=cmap,
                          extend=extend, norm=LogNorm(vmin=vmin, vmax=vmax))
    else:
        cf = plt.contourf(x, y, z, vmin=vmin, vmax=vmax, levels=levels,
                          cmap=cmap, extend=extend)

    # Station markers for reference
    plt.scatter(p2s.lons, p2s.lats, marker="v", alpha=0.75, c="None", 
                ec="w", s=15, zorder=8)

    # Annotate station name 
    if False:
        for lon_, lat_, id_ in zip(p2s.lons, p2s.lats, p2s.ids):
            plt.text(lon_, lat_, id_, fontsize=7, color="w", alpha=0.5)
        
    # Plot source as location or mechanism
    mt = [p2s.cmt["Mrr"], p2s.cmt["Mtt"], p2s.cmt["Mpp"], 
          p2s.cmt["Mrt"], p2s.cmt["Mrp"], p2s.cmt["Mtp"]]
    # Kludge past the fact that ObsPy can't plot the explosion-type nonisotropic
    # source correctly so we replace with an isotropic source
    if p2s.tag.endswith("013") or p2s.tag.endswith("NK6"):
        mt = [3.258319e+23, 3.25825e+23, 3.25818e+23, 
              24315570.0, -1093.535, -301.5596]
    bb = beach(mt, xy=(p2s.srclon, p2s.srclat), width=5e-2, 
               facecolor=mt_color, edgecolor="k", zorder=11, linewidth=1)
    ax.add_collection(bb)

    # Add coastline to the map (NK specific)
    if os.path.exists(fid_coastline):
        df = pd.read_csv(fid_coastline)
        for i in df["polygon_id"].unique():
            df[df["polygon_id"] == i].plot(x="longitude", y="latitude", 
                                           ax=ax, legend=False, c=coast_color, 
                                           lw=1.5, zorder=9
                                           )

    # Ruler(s): line + great-circle distance annotation between point pairs
    if ruler:
        # Accept either a single (point1, point2) pair or a list of pairs
        if not isinstance(ruler[0][0], (list, tuple, np.ndarray)):
            ruler = [ruler]
        for point1, point2 in ruler:
            add_ruler(ax, point1, point2)

    # Accoutremont
    plt.xlim([128.1, 129.78])  # Manually set
    plt.ylim([40.56, 42.325])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.25))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.25))
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.set_aspect("equal")

    plt.xlabel("Longitude", fontsize=14)
    plt.ylabel("Latitude", fontsize=14)
    if not title:
        tag, src = p2s.tag.split("_")
        src = src.split("-")[-1]
        tag = {"ALPHA": "1D", "BETA": "1D_TOPO", "CHARLIE": "1D_TOPO_SCTR"}[tag]

        title = (f"Model={tag}; Src={src}; Z={p2s.cmt['depth']}km; "
                 f"Comp={p2s.components}; Freq={p2s.fmin}\u2013{p2s.fmax} Hz")
    ax.set_title(title, fontsize=16)

    # Reserve a fixed-width strip on the right for the colorbar via `rect`,
    # rather than just letting tight_layout() use the full figure width.
    # This matters because the colorbar's own Axes (cax) is appended AFTER
    # tight_layout() runs (see below) -- if tight_layout() weren't told to
    # leave this margin, it would let `ax` stretch to fill the whole figure,
    # leaving no room for the colorbar's tick labels once cax gets carved
    # out of that space, clipping them off the right edge of the figure.
    # The reserved fraction is constant regardless of `colorbar`, so `ax`
    # still doesn't move when the colorbar is toggled on/off.
    plt.tight_layout(rect=[0, 0, 0.90, 1])

    # Colorbar lives on its own Axes, appended here -- AFTER tight_layout()
    # rather than before. tight_layout() recomputes subplot spacing from
    # each Axes' actual rendered content, so if the (empty-when-toggled-off)
    # colorbar Axes already existed, tight_layout() would size `ax` based on
    # whether a colorbar happened to be drawn in it. Appending it afterward
    # means `ax` is laid out once, from its own content only, and stays the
    # same size/position regardless of the `colorbar` toggle.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.125)
    if colorbar:
        cb = f.colorbar(cf, cax=cax, extend=extend, ticks=ticks)
        cb.set_label(rf"$\leftarrow$ Earthquake    "
                     rf"[P/S Ratio]     "
                     rf"Explosion $\rightarrow$",
                     fontsize=15, c="k")
        cb.ax.get_yaxis().labelpad = -50
        cb.ax.tick_params(labelsize=14)
        if log_scale:
            # Default LogNorm tick labels render in scientific notation
            # (e.g. "8x10^-1") -- show the same plain decimal values used
            # on a linear scale instead
            cb.ax.yaxis.set_minor_formatter(ticker.NullFormatter())
            cb.set_ticklabels([f"{t:g}" for t in ticks])
    else:
        cax.axis("off")

    # Subplot label -- `ax` is now the same size/position either way, so
    # this no longer needs a separate position per colorbar on/off
    if subplot_label:
        plt.figtext(0.09, .9, s=subplot_label, c="k", fontsize=27.5,
                    zorder=250, ha="left", va="top")

    # overlay waveforms + p/s picks for user-requested stations directly on
    # the map, pinned to their true lon/lat location. this must happen after
    # set_aspect()/tight_layout() above -- both can still move/resize `ax`,
    # and the inset boxes are placed in absolute figure-fraction coordinates
    # derived from ax.transdata, so computing them any earlier leaves the
    # insets stranded from their station markers once the layout settles.
    if highlight_stations:
        f.canvas.draw()  # finalize the layout so transdata is up to date
        for station in highlight_stations:
            result = p2s.get_station_waveform(station)
            if result is None:
                continue
            tr, p_idx, s_idx, p_window, s_window = result

            ids = list(p2s.ids)
            if station not in ids:
                print(f"station '{station}' not in processed station list, "
                      f"skipping map overlay")
                continue

            idx = ids.index(station)
            lon_, lat_ = p2s.lons[idx], p2s.lats[idx]

            add_waveform_inset(f, ax, lon_, lat_, tr, p_idx, s_idx,
                               width=.25, height=0.15, annotate=annotate)
            ax.scatter(lon_, lat_, marker="v", c="w", ec="k",  zorder=100,
                       s=30)

    if save:
        plt.savefig(save)
    if show:
        plt.show()


def parse_args():
    parser = argparse.ArgumentParser(description="Synthetic P/S Heatmaps")
    parser.add_argument("path", type=str, nargs="?", default=None,
                        help="path to waveform files .sem?")
    parser.add_argument("-C", "--choice", type=str, nargs="?", default="taup",
                        help="choice for windowing style")
    parser.add_argument("-c", "--components", default="ZNE", type=str, 
                        nargs="?", help="components to include")
    parser.add_argument("-f1", "--fmin", type=float, default=2, nargs="?",
                        help="minimum filter frequency")
    parser.add_argument("-f2", "--fmax", type=float, default=8, nargs="?",
                        help="maximum filter frequency")
    parser.add_argument("-p", "--parallel", action="store_true", 
                        help="use multiprocessing to evaluate in parallel")
    parser.add_argument("-n", "--ntasks", default=os.cpu_count(), nargs="?", 
                        type=int, help="how many cores to use")
    parser.add_argument("-o", "--overwrite", default=False, action="store_true",
                        help="overwrite existing data file")
    parser.add_argument("--colorbaroff", action="store_true")
    parser.add_argument("--log-scale", action="store_true",
                        help="color the heatmap on a log scale instead of "
                             "linear")
    parser.add_argument("--annotate", action="store_true")
    parser.add_argument("--save", type=str, nargs="?", default=None)
    parser.add_argument("--noshow", action="store_true")

    parser.add_argument("-i", default=0, type=int, nargs="?",
                        help="starting index for processing")
    parser.add_argument("-j", default=-1, type=int, nargs="?",
                        help="ending index for processing")

    parser.add_argument("-S", "--stations", type=str, nargs="+", default=None,
                        help="station codes to highlight on the heatmap; "
                             "each gets a marker plus an inset plot of its "
                             "waveform and P/S picks at its map location")
    parser.add_argument("--grid-res", type=int, default=200,
                        help="interpolation grid resolution (points per "
                             "axis) for the heatmap")
    parser.add_argument("--smooth-sigma", type=float, default=1.5,
                        help="Gaussian smoothing sigma (in grid cells) "
                             "applied to the interpolated heatmap; 0 to "
                             "disable")
    parser.add_argument("--ruler", type=float, nargs="+", default=None,
                        help="draw a ruler line + distance annotation (km) "
                             "between two points on the heatmap; give one "
                             "or more groups of 'lon1 lat1 lon2 lat2'")

    return parser.parse_args()


def main():
    args = parse_args()
   
    # Input directories/files
    path = Path(args.path)  # e.g. ALPHA/TDL-015
    tag = path.name
    model = path.parent

    # Custom set the subplot labels
    try:
        subplot_label = {"ALPHA/TDL-015": "A)",
                         "ALPHA/TDL-001": "B)",
                         "ALPHA/TDL-013": "C)",
                         "ALPHA/NK6": "D)",
                         "CHARLIE/TDL-015": "A)",
                         "CHARLIE/TDL-001": "B)",
                         "CHARLIE/TDL-013": "C)",
                         "CHARLIE/NK6": "D)",
                         }[str(path)]
    except KeyError:
        subplot_label = ""


    path_src = f"tour_de_lune/CMTSOLUTION_{tag}"
    path_sta = "STATIONS"

    # Output directories/files
    path_save = f"data/{model}_{tag}_{args.fmin}-{args.fmax}_{args.components}.npz"
    path_fig  = f"figures/{model}_{tag}/"

    # Run ratio maker
    p2s = P2SRatio(path=path, path_src=path_src, path_sta=path_sta,
                   path_save=path_save, path_fig=path_fig,
                   choice=args.choice,
                   fmin=args.fmin, fmax=args.fmax,
                   components=args.components, 
                   overwrite=args.overwrite,)
    p2s.calculate_ratio(i=args.i, j=args.j, parallel=args.parallel,
                        ntasks=args.ntasks)

    # Create rulers to add to heatmap figure
    ruler = None
    if args.ruler:
        if len(args.ruler) % 4 != 0:
            raise ValueError("--ruler must be given in groups of 4 floats: "
                             "lon1 lat1 lon2 lat2")
        ruler = [((args.ruler[i], args.ruler[i + 1]),
                  (args.ruler[i + 2], args.ruler[i + 3]))
                 for i in range(0, len(args.ruler), 4)]

    if args.save:
        save = args.save
    else:
        save = os.path.join("figures", f"heatmap_{model}_{tag}.png")

    plot_heatmap(p2s, save=save, show=not args.noshow,
                 method="linear", grid_res=args.grid_res,
                 smooth_sigma=args.smooth_sigma,
                 highlight_stations=args.stations, ruler=ruler,
                 subplot_label=subplot_label, colorbar=not args.colorbaroff,
                 annotate=args.annotate, log_scale=args.log_scale
                 )


if __name__ == "__main__": 
    main()
