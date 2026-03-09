"""
Calculate and scatter plot P-to-S amplitude ratio for a set of synthetics.
Timing windows based on TauP. Allows inputting multiple sets of synthetics
to compare different events etc.
"""
import os
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd


from pathlib import Path
from scipy.interpolate import griddata
from obspy.imaging.beachball import beach
from obspy.taup import TauPyModel
from obspy.geodetics import kilometers2degrees
from pysep import read_sem, logger

logger.setLevel("CRITICAL")

# For TauP arrival windows
P_PHASE_LIST = ["p", "P", "PP", "pP", "Pn", "Pg"]
S_PHASE_LIST = ["s", "S", "SS", "sS", "Sn", "Sg"]
TAUP_MODEL = "iasp91"

# For beachball plotting
MOMENT_TENSORS = {
    "NK6ISO": [5.395, 5.395, 5.395, 0, 0, 0],
    "EQ6": [3.54, 24, -27.59, 11.22, -4.015, 21.55],
    "NK1": [6.06, 3.18, 4.94, 1.32, 3.87, .375],
    "NK6": [7.14, 5.17, 4.27, .130, -2.59, .451]
}

# For titles
DESCRIPTORS = {
    "NK6ISO": "Isotropic Explosion @ 1km",
    "EQ6": "Tectonic Earthquake (So. Korea) @ 5km",
    "NK1": "Alvizuri & Tape (2018) MT #1 @ 1km",
    "NK6": "Alvizuri & Tape (2018) MT #6 @ 1km",
    "NK6_5km": "Alvizuri & Tape (2018) MT #6 @ 5km",

}

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


def get_taup_arrivals(source_depth_in_km, distance_in_km, model=TAUP_MODEL):
    """Get arrival time windows from TauP"""
    dist_deg = kilometers2degrees(distance_in_km)

    model = TauPyModel(model=model)
    p_arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km,
                                        distance_in_degree=dist_deg,
                                        phase_list=P_PHASE_LIST)
    if not p_arrivals:
        return None, None
    
    p_arrivals = [_.time for _ in p_arrivals]
    p_window = [min(p_arrivals), max(p_arrivals)]

    s_arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km,
                                        distance_in_degree=dist_deg,
                                        phase_list=S_PHASE_LIST)
    if not s_arrivals:
        return None, None
    
    s_arrivals = [_.time for _ in s_arrivals]
    s_window = [min(s_arrivals), max(s_arrivals)]

    return p_window, s_window


def plot_tr(tr, p_idx=None, s_idx=None, p_window=None, s_window=None,
            save="./", show=False):
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
    plt.title(f"{tr.get_id()}; dist={tr.stats.sac["dist"]:.2f}km; "
              f"baz={tr.stats.sac["baz"]:.2f}deg")

    if show:
        plt.show()
    if save:
        plt.savefig(os.path.join(save, f"{tr.get_id()}.png"))
    plt.close("all")


class P2SRatio:
    """Gather P2S ratio in class attributes"""
    def __init__(self, path, source="CMTSOLUTION", stations="STATIONS",
                 components="Z", plot_waveforms="./figures", overwrite=False):
        """Setup attributes for"""
        self.path = Path(path)
        self.tag = self.path.name
        self.path_save = self.path / f"data_{components}.npz"

        self.source = source
        self.stations = stations
        self.components = components
        self.overwrite = overwrite

        # Only use the first component to get IDs, even if we have multi-comp.
        self.fids = sorted(self.path.glob(f"*HX[{components[0]}].sem?"))

        if plot_waveforms:
            self.plot_waveforms = os.path.join(plot_waveforms, self.tag)
            os.makedirs(self.plot_waveforms, exist_ok=True)
        else:
            self.plot_waveforms = None

        self.p2sratios = []
        self.lats = []
        self.lons = []
        self.distances = []
        self.ids = []
        self.srclat = None
        self.srclon = None

    def __len__(self):
        return len(self.p2sratios)
    
    def write(self):        
        dict_out = {"p2sratios": np.array(self.p2sratios),
                    "lats": np.array(self.lats),
                    "lons": np.array(self.lons),
                    "distances": np.array(self.distances),
                    "ids": np.array(self.ids),
                    }
        np.savez(self.path_save, **dict_out)
    
    def read(self):
        dict_in = np.load(self.path_save)
        self.p2sratios = dict_in["p2sratios"]
        self.lats = dict_in["lats"]
        self.lons = dict_in["lons"]
        self.distances = dict_in["distances"]    
        # self.ids = dict_in["stations"]
    
    def get_data(self, fmin=2, fmax=4, i=0, j=-1):
        """Read or fetch data from files"""
        if self.path_save.exists() and (not self.overwrite):
            self.read()
        else:
            error = 0
            for fid in self.fids[i:j]:
                try:
                    self.get_data_single(fid, fmin, fmax)
                except Exception as e:
                    breakpoint()
                    error += 1
                    continue
            if error:
                print(f"{self.tag}: {error}/{len(self.fids[i:j])} errors")
                print(e)
            self.write()

        tr = read_sem(self.fids[0], source=self.source, 
                      stations=self.stations)[0]
        self.srclat = tr.stats.sac["evla"]
        self.srclon = tr.stats.sac["evlo"]

    def get_data_single(self, fid, fmin, fmax):
        """Read files and collect information including P2S"""
        tr = read_sem(fid, source=self.source, stations=self.stations)[0]
        tr.filter("bandpass", freqmin=fmin, freqmax=fmax)

        p_window, s_window = get_taup_arrivals(tr.stats.sac["evdp"], 
                                               tr.stats.sac["dist"])
        if not p_window or not s_window:
            if self.plot_waveforms:
                self.plot_tr(tr)
            raise Exception
        
        p2sratio, p_idx, s_idx = get_p2s(tr=tr, p_window=p_window, 
                                         s_window=s_window)  
        if self.plot_waveforms:
            plot_tr(tr, p_idx, s_idx, p_window, s_window, 
                    save=self.plot_waveforms)


        # If more than 1 component we need to get the other 2 arrays
        if len(self.components) > 1:
            comp = self.components[0]
            for comp_ in self.components[1:]:
                fid_ = str(fid).replace(f"X{comp}.sem", f"X{comp_}.sem")
                tr_ = read_sem(fid_, source=self.source, 
                               stations=self.stations)[0]
                tr_.filter("bandpass", freqmin=fmin, freqmax=fmax)

                p2sratio_, p_idx_, s_idx_ = get_p2s(tr=tr_, p_window=p_window, 
                                                    s_window=s_window)  
                p2sratio += p2sratio_
                if self.plot_waveforms:
                    plot_tr(tr_, p_idx_, s_idx_, p_window, s_window,
                            save=self.plot_waveforms)
            
            p2sratio /= len(self.components)  # take average

        self.p2sratios.append(p2sratio)
        self.ids.append(tr.stats.station)
        self.lats.append(tr.stats.sac["stla"])
        self.lons.append(tr.stats.sac["stlo"])


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


def plot_heatmap(p2s, fmin=1, fmax=6, threshold=0.8):
    """
    For a single event plot a map of amplitude ratios for each station 
    interpolated to create a continuous figure rather than a scatterplot.
    """
    # Custom parameters here
    cmap = "seismic"
    mt_color = "gold"
    coast_color = "k"
    fid_coastline = "coastline_128_130_40_43.csv"
    # ^^^

    f, ax = plt.subplots(1, figsize=(15, 12))

    # Create a regular grid from lon/lat ranges
    lon_grid = np.linspace(min(p2s.lons), max(p2s.lons), 
                            2*len(np.unique(p2s.lons)))
    lat_grid = np.linspace(min(p2s.lats), max(p2s.lats), 
                            2*len(np.unique(p2s.lats)))
    x, y = np.meshgrid(lon_grid, lat_grid)
    
    # Interpolate p2s ratios onto the grid
    # Kludge to force all ratios to be positive so that we don't get negative
    # values. They should be positive but for some reason some are not
    # z_vals = np.abs(p2s.p2sratios)
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
    plt.scatter(p2s.lons, p2s.lats, marker="v", alpha=0.25, c="None", 
                ec="w", s=10, zorder=8)
        
    # Plot source as location or mechanism
    mt = MOMENT_TENSORS[p2s.tag.split("_")[0]]
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
    plt.title(f"{DESCRIPTORS[p2s.tag]}; comp={p2s.components}; "
              f"freq={fmin}-{fmax} Hz; Threshold={threshold}", 
              fontsize=16)

    plt.tight_layout()
    plt.savefig(f"./figures/{p2s.tag}_heatmap.png")
    plt.show()


def main(path, components, fmin, fmax, overwrite, i=0, j=-1):
    """Get data and plot"""
    p2s = P2SRatio(path, source=f"CMTSOLUTION_{Path(path).name}", 
                components=components, overwrite=overwrite)
    p2s.get_data(fmin=fmin, fmax=fmax, i=i, j=j)
    plot_heatmap(p2s, fmin, fmax)


if __name__ == "__main__": 
    main(path=sys.argv[1], components="ZNE", fmin=1, fmax=6, overwrite=False)

