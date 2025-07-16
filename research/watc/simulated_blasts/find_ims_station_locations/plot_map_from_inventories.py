"""
Basic Cartopy plotter for ObsPy Inventory objects
"""
import os
import argparse
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from obspy import read_inventory


# Used for coordinate transforms
REF_PROJ = ccrs.PlateCarree()

def parse_args():
    """
    arg parsarg arg arg!
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files", type=str, nargs="+",
                        help="Inventory files to plot")
    parser.add_argument("-s", "--highlight", type=str, nargs="+",
                        help="Plot certain stations with different traits")
    parser.add_argument("-d", "--dpi", type=float, default=300)
    parser.add_argument("-l", "--event_latlon", type=tuple, 
                        default=(43.3430, 129.0360), 
                        help="Put an event location, fmt = (network.station)")
    parser.add_argument("-o", "--output", type=str, default="figures/map.png")

    return parser.parse_args()


def setup(proj_str="Stereographic", central_longitude=-147, central_latitude=0,
          figsize=None, dpi=None, lw_axis=1.5, lw_coast=0.5,
          default_proj_str="Stereographic"):
    """
    Setup a standard looking Cartopy global plot with coastline and gridlines

    :type proj_str: str
    :param proj_str: The string representing the Cartopy CRS projection 
        (default: "Stereographic")
    :type central_longitude: float
    :param central_longitude: The central longitude of the projection 
        (default: 0)
    :type central_latitude: float
    :param central_latitude: The central latitude of the projection 
        (default: 0)
    :type figsize: tuple
    :param figsize: The size of the figure (default: None)
    :type dpi: int
    :param dpi: The resolution of the figure in dots per inch (default: None)
    :type axis_linewidth: float
    :param axis_linewidth: The linewidth of the axes (default: 1.5)
    :type default_proj_str: str
    :param default_proj_str: The default string representing the Cartopy CRS 
        projection (default: "Stereographic")
    :rtype: tuple
    :return: The created figure and axes objects
    """
    # Attempt to grab a plate projection from Cartopy
    try:
        projection = getattr(ccrs, proj_str)
    except AttributeError:
        print(f"{proj_str} is not a valid Cartopy CRS projection, setting "
              f"default value: '{default_proj_str}'")
        projection = getattr(ccrs, default_proj_str)

    # Most major projections use the central_longitude function while only
    # some require central latitude, so we set that outside the init
    projection = projection(central_longitude=central_longitude)
    projection.central_latitude = 0

    # Start fig either standalone or tandem with a waveform plot in GridSpec
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = plt.axes(projection=projection)

    ax.coastlines(lw=lw_coast, color="darkgray", zorder=1)
    ax.add_feature(cartopy.feature.OCEAN, zorder=0)
    ax.add_feature(cartopy.feature.LAND, zorder=0)

    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False,
                      y_inline=False, linewidth=.25, alpha=0.25, color="k")

    gl.top_labels = True
    gl.right_labels = False
    gl.left_labels = True
    gl.bottom_labels = False
    gl.xlabel_style = {"fontsize": 5, "rotation": 0}
    gl.ylabel_style = {"fontsize": 5, "rotation": 0}

    # Axis linewidth is set differently than in Matplotlib, see:
    ax.spines["geo"].set_linewidth(lw_axis)

    # Remove whitespace
    plt.tight_layout()

    return fig, ax


def plot_inventory(inv, ax, color="r", marker="v", size=10., alpha=1, 
                   fontsize=2.75, fontcolor="k", text=True, zorder=1, 
                   highlight=None, plotted=None):
    """
    Plot the stations in an ObsPy Inventory object
    """
    # We will use this to make sure we don't double plot stations
    if not plotted:
        plotted = []

    for network in inv:
        for station in network:
            lon = station.longitude
            lat = station.latitude
            code = f"{network.code}.{station.code}"

            # Highlight requested stations
            if highlight and code in highlight:
                ax.scatter(lon, lat, color=color, marker=marker, ec="k", 
                           lw=1.1, s=size + 2, alpha=alpha, transform=REF_PROJ, 
                           zorder=zorder)
                plotted.append(code)

            if code not in plotted:
                ax.scatter(lon, lat, color=color, marker=marker, 
                           s=size, alpha=alpha, transform=REF_PROJ, 
                           zorder=zorder, ec="k", lw=0.25)
                plotted.append(code)

            if text:
                # Right edge of the map, keep text in the map
                if lon >= 165. and lon < 180.:
                    horizontalalignment = "right"
                else:
                    horizontalalignment = "left"
                ax.text(lon, lat, s=f"   {network.code}.{station.code}", 
                        transform=REF_PROJ, fontsize=fontsize, c=fontcolor, 
                        # bbox=dict(facecolor="w", edgecolor="k", pad=0.5),
                        horizontalalignment=horizontalalignment,
                        zorder=zorder + 1)
    return ax, plotted


def plot_event(ev_lat, ev_lon, ax, color="y", marker="*", size=40, alpha=1):
    """
    Plot an event as a yellow star or whatever
    """
    ax.scatter(ev_lon, ev_lat, ec="k", lw=0.5, color=color, marker=marker, 
               s=size, alpha=alpha, transform=REF_PROJ)


def plot_test_sites(ax):
    """
    Plot some test site locations on the map
    """
    sites = {
            "Lop Nur": (40.778373, 89.264008),
            "NK": (41.343, 129.036),
            "Novaya Zemlya": (73.4, 54.9),
            "Semey": (50.433333, 80.266667),
            "Moruroa": (-21.833333, -138.833333),
            "Pokhran": (27.078889, 71.722500),
            "Mushaf": (32.048611, 72.665278),
            "Ras Koh": (28.828531, 65.194953),
            "NNSS": (37.116667, -116.050000),
            "Amchitka": (51.542222, 178.983333),
            "Kiritimati": (1.788280, -157.403378),
            "Emu": (-29.88980, 131.64670),
            "In Ekker": (4.74000, 23.93000),

            }
    for name, loc in sites.items():
        lat, lon = loc
        ax.scatter(lon, lat, color="yellow", marker="*", ec="k",
                   lw=.25, s=13, transform=REF_PROJ, zorder=10)
        # Right edge of the map, keep text in the map
        if lon >= 165. and lon < 180.:
            horizontalalignment = "right"
        else:
            horizontalalignment = "left"
        ax.text(lon, lat, s=f"  {name}", transform=REF_PROJ, fontsize=3,
                c="k", zorder=11, horizontalalignment=horizontalalignment
                )



if __name__ == "__main__":
    args = parse_args()
    ev_lat, ev_lon = args.event_latlon

    # Setting kwargs
    inv_dict = {}
    for i, f in enumerate(args.files):
        inv_dict[f] = {"inv": read_inventory(f),
                       "color": f"C{i}",
                       "marker": "s",
                       "size": 10.,
                       "fontsize": 2.75,
                       "fontcolor": "k",
                       }
        if "primary" in f:
            inv_dict[f]["color"] = "darkviolet"
        elif "auxiliary" in f:
            inv_dict[f]["color"] = "dodgerblue"
        elif "WATC" in f:
            inv_dict[f]["color"] = "red"
            inv_dict[f]["marker"] = "v"
            inv_dict[f]["zorder"] = 10
        elif "PRIMARY" in f:
            inv_dict[f]["color"] = "magenta"
            inv_dict[f]["marker"] = "^"
            inv_dict[f]["size"] = 10
            # inv_dict[f]["text"] = False
        elif "AUXILIARY" in f:
            inv_dict[f]["color"] = "powderblue"
            inv_dict[f]["marker"] = "^"
            inv_dict[f]["size"] = 10
            # inv_dict[f]["text"] = False

    # Plotting
    plotted = []
    fig, ax = setup(proj_str="Robinson", dpi=args.dpi)
    for key, vals in inv_dict.items():
        inv = vals.pop("inv")
        ax, plotted = plot_inventory(inv=inv, ax=ax, highlight=args.highlight, 
                                     plotted=plotted, **vals)

    # plot_event(ev_lat=ev_lat, ev_lon=ev_lon, ax=ax)
    plot_test_sites(ax)

    # Finalize
    if not os.path.exists(os.path.dirname(args.output)):
        os.mkdir(os.path.dirname(args.output))
    plt.savefig(args.output, dpi=args.dpi)
    plt.show()
