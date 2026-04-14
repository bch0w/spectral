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
from scipy.interpolate import griddata
from obspy.imaging.beachball import beach
from obspy.taup import TauPyModel
from obspy.geodetics import kilometers2degrees
from pysep import read_sem, logger

logger.setLevel("CRITICAL")


FID_COAST = "/home/bhchow/work/data/cartography/coastline_128_130_40_43.csv"

def get_p2s(tr, p_window, s_window, choice="before_after_s"):
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
    s_start = int(s_start)
    s_stop = int(s_stop)

    # Only look inside the given P or S window
    if choice == "window":    
        p_idx = np.argmax(tr[p_start: p_stop])
        s_idx = np.argmax(tr[s_start: s_stop]) + s_start
    # Use the extremum bounds to search for P and S
    elif choice == "bounds":
        p_idx = np.argmax(tr[:p_stop])
        s_idx = np.argmax(tr[s_start:]) + s_start
    # Use min S time as a dividing line
    elif choice == "before_after_s":
        p_idx = np.argmax(tr[:s_start])
        s_idx = np.argmax(tr[s_start:]) + s_start 

    # P-to-S amplitude ratio, index of max P, index of max S
    return float(np.abs(tr[p_idx] / tr[s_idx])), p_idx, s_idx


def get_taup_arrivals(source_depth_in_km, distance_in_km, p_phase_list=None,
                      s_phase_list=None, model="iasp91"):
    """
    Get arrival time windows from TauP for a given `TAUP_MODEL`
    """
    # By default query all the crustal directarrivals
    if not p_phase_list:
        p_phase_list = ["p", "P", "PP", "pP", "Pn", "Pg"]
    if not s_phase_list:
        s_phase_list = ["s", "S", "SS", "sS", "Sn", "Sg"]

    dist_deg = kilometers2degrees(distance_in_km)

    model = TauPyModel(model=model)
    p_arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km,
                                        distance_in_degree=dist_deg,
                                        phase_list=p_phase_list)
    if not p_arrivals:
        return None, None
    
    p_arrivals = [_.time for _ in p_arrivals]
    p_window = [min(p_arrivals), max(p_arrivals)]

    s_arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km,
                                        distance_in_degree=dist_deg,
                                        phase_list=s_phase_list)
    if not s_arrivals:
        return None, None
    
    s_arrivals = [_.time for _ in s_arrivals]
    s_window = [min(s_arrivals), max(s_arrivals)]

    return p_window, s_window


def plot_tr(tr, p_idx=None, s_idx=None, p_window=None, s_window=None,
            title_prepend="", save="./", show=False):
    """Make individual waveform plot"""
    f, ax = plt.subplots(1, figsize=(8,6), dpi=100)
    plt.plot(tr.times(), tr.data, c="k", zorder=9, lw=1)

    if p_idx:
        plt.scatter(tr.times()[p_idx], tr.data[p_idx], zorder=10, 
                    c="r", label=f"P {tr.data[p_idx]:.2E}")
    if s_idx:
        plt.scatter(tr.times()[s_idx], tr.data[s_idx], zorder=10,
                    c="g", label=f"S {tr.data[s_idx]:.2E}")
    if p_window:
        plt.axvline(p_window[0], c="r", zorder=8, lw=1)
        plt.axvline(p_window[1], c="r", zorder=8, lw=1)
    if s_window:
        plt.axvline(s_window[0], c="g", zorder=8, lw=1)
        plt.axvline(s_window[1], c="g", zorder=8, lw=1)

    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel("Velocity [m/s]")
    plt.title(f"{title_prepend}"
              f"{tr.get_id()}; dist={tr.stats.sac["dist"]:.2f}km; "
              f"baz={tr.stats.sac["baz"]:.2f}deg"
              )

    if show:
        plt.show()
    if save:
        plt.savefig(os.path.join(save, f"{tr.get_id()}.png"))
    plt.close(f)


class P2SRatio:
    """Gather P2S ratio in class attributes"""
    def __init__(self, path, path_src="CMTSOLUTION", path_sta="STATIONS", 
                 path_save="data.npz", path_fig="figures", fmin=1, fmax=6, 
                 components="ZNE", overwrite=False):
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
        self.path = Path(path)
        self.tag = self.path.name

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
        self.fids = sorted(self.path.glob(f"*X[{components[0]}].sem?"))
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
        return read_sem(fid, source=self.path_src, stations=self.path_sta)[0]

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
            self.read()
        else:
            if not parallel:
                print(f"calculating p2s ratios in serial for {len(self.fids)}")
                error = 0
                for fid in self.fids[i:j]:
                    try:
                        tr, p2sratio = self.calculate_ratio_single(fid)
                        self.p2sratios.append(p2sratio)
                        self.ids.append(tr.stats.station)
                        self.lats.append(tr.stats.sac["stla"])
                        self.lons.append(tr.stats.sac["stlo"])
                    except Exception as e:
                        error += 1
                        continue
                if error:
                    print(f"{self.tag}: {error}/{len(self.fids[i:j])} errors")
                    print(e)
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
                    try:
                        tr, p2sratio = future.result()
                        self.p2sratios.append(p2sratio)
                        self.ids.append(tr.stats.station)
                        self.lats.append(tr.stats.sac["stla"])
                        self.lons.append(tr.stats.sac["stlo"])
                    except Exception as e:
                        continue

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
        tr = self._read_sem(fid)
        tr.filter("bandpass", freqmin=self.fmin, freqmax=self.fmax)

        p_window, s_window = get_taup_arrivals(tr.stats.sac["evdp"], 
                                               tr.stats.sac["dist"])
        if not p_window or not s_window:
            if self.path_fig:
                self.plot_tr(tr)
            raise Exception
        
        p2sratio, p_idx, s_idx = get_p2s(tr=tr, p_window=p_window, 
                                         s_window=s_window)  
        if self.path_fig:
            plot_tr(tr, p_idx, s_idx, p_window, s_window, 
                    title_prepend=f"{self.tag} ",
                    save=self.path_fig)

        # If more than 1 component we need to get the other 2 arrays
        if len(self.components) > 1:
            comp = self.components[0]
            for comp_ in self.components[1:]:
                fid_ = str(fid).replace(f"X{comp}.sem", f"X{comp_}.sem")
                tr_ = self._read_sem(fid_)
                tr_.filter("bandpass", freqmin=self.fmin, freqmax=self.fmax)

                p2sratio_, p_idx_, s_idx_ = get_p2s(tr=tr_, p_window=p_window, 
                                                    s_window=s_window)  
                p2sratio += p2sratio_
                if self.path_fig:
                    plot_tr(tr_, p_idx_, s_idx_, p_window, s_window,
                            title_prepend=f"{self.tag} ",
                            save=self.path_fig)
            
            p2sratio /= len(self.components)  # take average

        return tr, p2sratio


def plot_scatterplot(paths, fmin=2, fmax=4, components="Z", j=-1):
    """
    Run through P2S for multiple events and create a scatterplot of distance
    versus P2S amplitude ratio. Similar to Walter et al 2018 F2a but without
    magnitude on the x axis
    """
    ratios = []
    for path in paths:
        p2s = P2SRatio(path, source=f"CMTSOLUTION_{Path(path).name}", 
                       components=components)
        print(p2s.tag)
        p2s.get_data(fmin=fmin, fmax=fmax, j=j)
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


def plot_heatmap(p2s, threshold=0.8, save="./figures", cmap="seismic", 
                 mt_color="gold", coast_color="k", 
                 fid_coastline=FID_COAST):
    """
    For a single event plot a map of amplitude ratios for each station 
    interpolated to create a continuous figure rather than a scatterplot.
    """
    print("plotting heatmap")
    f, ax = plt.subplots(1, figsize=(15, 12))

    # Create a regular grid from lon/lat ranges
    lon_grid = np.linspace(min(p2s.lons), max(p2s.lons), 
                            2*len(np.unique(p2s.lons)))
    lat_grid = np.linspace(min(p2s.lats), max(p2s.lats), 
                            2*len(np.unique(p2s.lats)))
    x, y = np.meshgrid(lon_grid, lat_grid)
    
    # Interpolate p2s ratios onto the grid
    z = griddata((p2s.lons, p2s.lats), p2s.p2sratios, (x, y), 
                 method="linear")
    
    # Setting hard bounds for colorscale to keep all figures looking similar
    vmin = 0
    vmax = vmin + threshold * 2
    levels = np.linspace(vmin, vmax, int((vmax - vmin) * 10 * 6))
    ticks = [vmin, threshold, vmax]

    if vmax < z.max():
        extend = "max"
    else:
        extend = "neither"

    cf = plt.contourf(x, y, z, vmin=vmin, vmax=vmax, levels=levels, cmap=cmap,
                      extend=extend)
    
    # Colorbar that respects the vmin vmax values 
    cb = f.colorbar(cf, ax=ax, pad=0, extend=extend, ticks=ticks)
    cb.set_label(rf"$\leftarrow$ Earthquake    "
                 rf"[P/S Ratio]     "
                 rf"Explosion $\rightarrow$", 
                 fontsize=15, c="k")
    cb.ax.get_yaxis().labelpad = -60
    cb.ax.tick_params(labelsize=14)

    # Station markers for reference
    if False:
        plt.scatter(p2s.lons, p2s.lats, marker="v", alpha=0.25, c="None", 
                    ec="w", s=10, zorder=8)

        for lon_, lat_, id_ in zip(p2s.lons, p2s.lats, p2s.ids):
            plt.text(lon_, lat_, id_, fontsize=4, color="w")
        
    # Plot source as location or mechanism
    mt = [p2s.cmt["Mrr"], p2s.cmt["Mtt"], p2s.cmt["Mpp"], 
          p2s.cmt["Mrt"], p2s.cmt["Mrp"], p2s.cmt["Mtp"]]
    bb = beach(mt, xy=(p2s.srclon, p2s.srclat), width=5e-2, 
                facecolor=mt_color, edgecolor="k", zorder=11)
    ax.add_collection(bb)

    # Add coastline to the map (NK specific)
    if os.path.exists(fid_coastline):
        df = pd.read_csv(fid_coastline)
        for i in df["polygon_id"].unique():
            df[df["polygon_id"] == i].plot(x="longitude", y="latitude", 
                                           ax=ax, legend=False, c=coast_color, 
                                           lw=1.5, zorder=9
                                           )

    # Accoutremont
    plt.xlim([x.min(), x.max()])
    plt.ylim([y.min(), y.max()])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.25))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.25))
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.set_aspect("equal")

    plt.xlabel("Longitude (deg)", fontsize=14)
    plt.ylabel("Latitude (deg)", fontsize=14)
    plt.title(f"{p2s.tag} {p2s.cmt['depth']}km; "
              f"comp={p2s.components}; freq={p2s.fmin}-{p2s.fmax} Hz; "
              f"Threshold={threshold}", 
              fontsize=16)

    plt.tight_layout()
    if save:
        plt.savefig(save)
    plt.show()


def parse_args():
    parser = argparse.ArgumentParser(description="Synthetic P/S Heatmaps")
    parser.add_argument("-m", "--model", type=str, nargs="?", required=True,
                        # choices=["alpha", "beta", "charlie", "echo"],      
                        help="model options")
    parser.add_argument("-s", "--source", type=str, nargs="?",  required=True,
                        # choices=["EQ2", "NK6", "NK6b", "ISO"], 
                        help="source name")
    parser.add_argument("-c", "--components", default="ZNE", type=str, 
                        nargs="?", help="components to include")
    parser.add_argument("-f1", "--fmin", type=float, default=1, nargs="?",
                        help="minimum filter frequency")
    parser.add_argument("-f2", "--fmax", type=float, default=6, nargs="?",
                        help="maximum filter frequency")
    parser.add_argument("-p", "--parallel", action="store_true", 
                        help="use multiprocessing to evaluate in parallel")
    parser.add_argument("-n", "--ntasks", default=os.cpu_count(), nargs="?", 
                        type=int, help="how many cores to use")
    parser.add_argument("-o", "--overwrite", default=False, action="store_true",
                        help="overwrite existing data file")

    parser.add_argument("-i", default=0, type=int, nargs="?", 
                        help="starting index for processing")
    parser.add_argument("-j", default=-1, type=int, nargs="?", 
                        help="ending index for processing")

    return parser.parse_args()


def main():
    args = parse_args()
   
    # Input directories/files
    path = f"waveforms/{args.model}/{args.source}" 
    path_src  = (f"/home/bhchow/REPOS/spectral/research/watc/simblast/"
                 f"SPECFEM_DATA/CMTSOLUTIONS/paper_events/"
                 f"CMTSOLUTION_{args.source}")
    path_sta = (f"/home/bhchow/REPOS/spectral/research/watc/simblast/"
                f"SPECFEM_DATA/STATIONS/STATIONS_PAPER_NK_GRID")

    # Output directories/files
    path_save = (f"data/{args.model}_{args.source}_{args.fmin}_"
                 f"{args.fmax}_{args.components}.npz")
    path_fig  = f"figures/{args.model}/{args.source}"

    # Run ratio maker
    p2s = P2SRatio(path=path, path_src=path_src, path_sta=path_sta,
                   path_save=path_save, path_fig=path_fig,
                   fmin=args.fmin, fmax=args.fmax,
                   components=args.components, 
                   overwrite=args.overwrite,
                   )
    p2s.calculate_ratio(i=args.i, j=args.j, parallel=args.parallel,
                        ntasks=args.ntasks)

    plot_heatmap(p2s, save=f"figures/{args.model}_{args.source}.png")


if __name__ == "__main__": 
    main()
