"""
Make a nice pub quality map top be used for stuff
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from obspy.imaging.beachball import beach
from netCDF4 import Dataset
from matplotlib.ticker import MultipleLocator



def set_plot_aesthetic(
        ax, ytick_fontsize=12., xtick_fontsize=12., tick_linewidth=1.5,
        tick_length=5., tick_direction="in", ytick_format=None,
        xlabel_fontsize=12., ylabel_fontsize=12., axis_linewidth=1.5,
        spine_zorder=8, title_fontsize=14.,
        spine_top=True, spine_bot=True, spine_left=True, spine_right=True,
        xtick_minor=15, xtick_major=30, ytick_minor=15, ytick_major=30,
        xgrid_major=True, xgrid_minor=True, ygrid_major=True, ygrid_minor=True,
        **kwargs):
    """
    Set a uniform look for figures, stolen from PySEP
    """
    ax.title.set_fontsize(title_fontsize)
    ax.tick_params(axis="both", which="both", width=tick_linewidth,
                        direction=tick_direction, length=tick_length)
    ax.tick_params(axis="x", labelsize=xtick_fontsize)
    ax.tick_params(axis="y", labelsize=ytick_fontsize)
    ax.xaxis.label.set_size(xlabel_fontsize)
    ax.yaxis.label.set_size(ylabel_fontsize)

    # Thicken up the bounding axis lines
    for axis, flag in zip(["top", "bottom", "left", "right"],
                          [spine_top, spine_bot, spine_left, spine_right]):
        # Deal with the case where command line users are inputting strings
        if isinstance(flag, str):
            flag = bool(flag.capitalize() == "True")
        ax.spines[axis].set_visible(flag)
        ax.spines[axis].set_linewidth(axis_linewidth)

    # Set spines above azimuth bins
    for spine in ax.spines.values():
        spine.set_zorder(spine_zorder)

    # Scientific format for Y-axis
    if ytick_format == "sci":
        try:
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        except AttributeError:
            # If we are in log format axis this will not work
            pass
    else:
        try:
            ax.ticklabel_format(axis="y", style=ytick_format)
        except AttributeError:
            pass

    # Set xtick label major and minor which is assumed to be a time series
    if xtick_major:
        ax.xaxis.set_major_locator(MultipleLocator(float(xtick_major)))
    if xtick_minor:
        ax.xaxis.set_minor_locator(MultipleLocator(float(xtick_minor)))
    if ytick_minor:
        ax.yaxis.set_major_locator(MultipleLocator(float(ytick_major)))
    if ytick_major:
        ax.yaxis.set_minor_locator(MultipleLocator(float(ytick_minor)))

    plt.sca(ax)
    if xgrid_major:
        plt.grid(visible=True, which="major", axis="x", alpha=0.2, linewidth=1)
    if xgrid_minor:
        plt.grid(visible=True, which="minor", axis="x", alpha=0.2, linewidth=.5)
    if ygrid_major:
        plt.grid(visible=True, which="major", axis="y", alpha=0.2, linewidth=1)
    if ygrid_minor:
        plt.grid(visible=True, which="minor", axis="y", alpha=0.2, linewidth=.5)


def read_ncfile(fid):
    """
    Read in the ncfile and return its outputs
    """
    model  = Dataset(fid)
    latitude = model.variables["latitude"]
    longitude = model.variables["longitude"]
    rel_topo  = model.variables["relative_topo_radius"]

    return longitude, latitude, rel_topo


def plot_ncfile(fid, cmap="cividis", ncolor=43, vmin=-2, vmax=5, figsize=(16,9),
                dpi=150, anno_features=True, anno_landers=True, 
                anno_stations=None, anno_source=None, save=True, show=True):
    """
    Read in a target ncfile and plot the values in there for confirmation.
    Assuming parameter keys are the same as the written file
    """
    x, y, z = read_ncfile(fid)
    X, Y = np.meshgrid(x, y)

    f, ax = plt.subplots(1, figsize=figsize, dpi=dpi)

    cmap = plt.get_cmap(cmap, ncolor)
    pc = plt.pcolor(X, Y, z[:] * 1E-3, cmap=cmap, vmin=vmin, vmax=vmax)  # m -> km

    # Colorbar
    cbar = plt.colorbar(pc, fraction=0.024, pad=0.008, extend="both")
    cbar.set_label(label="Elevation [km]", size=12)
    cbar.outline.set_linewidth(1.5)

    set_plot_aesthetic((ax))

    # Text labels
    plt.title("VenusTopo719 Topography (R=6159km)")
    plt.xlabel("Longitude [$^{\circ}$]")
    plt.ylabel("Latitude [$^{\circ}$]")
    
    ax.set_aspect("equal")

    # Plot stations from STATIONS file
    if anno_stations:
        highlight = ["1310", "0909", "0817"]
        sta, net, lats, lons, d, e = np.loadtxt(anno_stations, dtype=str).T
        lats = [float(lat) for lat in lats]
        lons = [float(lon) for lon in lons]
        for i in range(len(sta)):
            if sta[i] in highlight:
                markersize = 5.5
                color = "forestgreen"
                alpha = 1.0
                fontsize = 8
                txt = f"  {sta[i]}"
            else:
                markersize = 2
                color = "w"
                alpha = 0.5
                fontsize = 5.5
                txt = sta[i]
            plt.plot(lons[i], lats[i], f"wv", markersize=markersize, 
                     markeredgecolor=color, lw=0.5, alpha=alpha)
            plt.text(lons[i], lats[i], txt,
                     # f"{net[i]}.{sta[i]}\n   {e[i]}",  # include elevation
                     c="w", 
                     fontsize=fontsize,
                     alpha=alpha)

    # Annotate source based on given latlon
    if anno_source:
        lat, lon = anno_source
        mt = [-1.56, -0.138, 1.69, -0.742, -0.121, -0.294]  # E26 dyn*cm
        bball = beach(fm=mt, xy=(lon, lat), facecolor="r", width=3.5, 
                      linewidth=0.5)
        ax.add_collection(bball)
        # plt.plot(lon, lat, "y*", markersize=10., markeredgecolor="k", lw=1)

    # Annotate features
    if anno_features:
        anno_features = {"Ishtar Terra":  (-25., 77.5),
                         "Beta Regio": (-105, 20.), 
                         "Alta Regio": (-162, -2.),  
                         "Alpha Regio": (3.6, -30.),
                         "Aphrodite Terra": (92.5, -20.),
                         }
        for name, coord in anno_features.items():
            plt.text(coord[0], coord[1], name, c="w", fontsize=13)# ,
                     #fontweight="bold")

    # Annotate lander locations, from Wikipedia
    # new lines are a hacky way to move the annotation without moving the marker
    if anno_landers:
        anno_landers = {"Venera 7\n": (-9, -5),
                        "Venera 8": (-25, -10),
                        "Venera 9": (-69, 31),
                        "Venera 10\n": (-69, 15.7),
                        "\nVenera 11": (-61, -14),
                        "Venera 12\n": (-66, -7),
                        "Venera 13\n": (-57, -7.08),
                        "\nVenera 14": (-50, -13.41),
                        "Vega 1": (177.11, 7.08),
                        "Vega 2": (177.11, -8.08),
                        }
        
        # alignments (ha, va)
        has = {"Vega 1": ("right", "bottom"),
               "Vega 2": ("right", "top"), 
               "\nVenera 11": ("right", "top"),
               "Venera 12\n": ("right", "bottom"),
               "\nVenera 14": ("left", "top"),
               }
        for name, coord in anno_landers.items():
            if name in has:
                ha, va = has[name]
            else:
                ha, va = ("left", "bottom")
            plt.plot(coord[0], coord[1], "wo", markersize=3, 
                     markeredgecolor="c", lw=1.0)
            plt.text(coord[0], coord[1], name, c="w", fontsize=8, ha=ha, va=va)

    f.tight_layout()

    if save:
        if isinstance(save, str):
            fid_out = save
        else:
            fid_out = os.path.basename(fid).split(".")[0] + ".png"
        plt.savefig(fid_out, transparent=True)

    if show:
        plt.show()


if __name__ == "__main__":
    choice = "grid"

    if choice == "map":
        plot_ncfile(fid="venus_topo_lat720_lon1440.nc",
                    cmap="cividis", ncolor=43, 
                    vmin=-2, vmax=5, 
                    figsize=(16,9), dpi=150, 
                    anno_features=True, 
                    anno_landers=True,
                    anno_stations="STATIONS",
                    anno_source=False,
                    show=True, save="venus_feature_map.png"
                    )
    elif choice == "grid":
        plot_ncfile(fid="venus_topo_lat720_lon1440.nc",
                    cmap="inferno", ncolor=23, 
                    vmin=-2, vmax=5, 
                    figsize=(16,9), dpi=200, 
                    anno_features=True, 
                    anno_landers=True,
                    anno_stations="STATIONS_GRID",
                    anno_source=(15, -78),  # (lat, lon)
                    show=True, save="venus_grid_map.png"
                    )





