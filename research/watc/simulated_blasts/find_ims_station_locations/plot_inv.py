"""
Basic Cartopy plotter for ObsPy Inventory objects
"""
import sys
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from obspy import read_inventory


# Used for coordinate transforms
REF_PROJ = ccrs.PlateCarree()


def setup(proj_str="Stereographic", central_longitude=0, central_latitude=0,
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

    ax.coastlines(lw=lw_coast, zorder=1)
    ax.add_feature(cartopy.feature.OCEAN, zorder=0)
    ax.add_feature(cartopy.feature.LAND, zorder=0)

    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False,
                      y_inline=False, linewidth=.25, alpha=0.25, color="k")

    gl.top_labels = False
    gl.right_labels = False
    # gl.left_labels = False
    # gl.bot_labels = False
    gl.xlabel_style = {"fontsize": 5, "rotation": 0}
    gl.ylabel_style = {"fontsize": 5, "rotation": 0}

    # Axis linewidth is set differently than in Matplotlib, see:
    ax.spines["geo"].set_linewidth(lw_axis)

    return fig, ax


def plot_inventory(inv, ax, color="r", marker="v", size=15, alpha=1, fontsize=5,
                   fontcolor="k", text=True):
    """
    Plot the stations in an ObsPy Inventory object
    """
    for network in inv:
        for station in network:
            lon = station.longitude
            lat = station.latitude
            ax.scatter(lon, lat, color=color, marker=marker, 
                       s=size, alpha=alpha, transform=REF_PROJ, zorder=10)
            if text:
                ax.text(lon, lat, s=f"   {network.code}.{station.code}", 
                        transform=REF_PROJ, fontsize=fontsize, c=fontcolor, 
                        # bbox=dict(facecolor="w", edgecolor="k", pad=0.5),
                        zorder=11)
    return ax


if __name__ == "__main__":
    files = sys.argv[1:]
    inv_dict = {}
    for i, f in enumerate(files):
        inv_dict[f] = {"inv": read_inventory(f),
                       "color": f"C{i}",
                       "marker": "s",
                       "size": 12.5,
                       "fontsize": 2.75,
                       "fontcolor": "k"
                       }
        if "primary" in f:
            inv_dict[f]["color"] = "darkviolet"
        elif "auxiliary" in f:
            inv_dict[f]["color"] = "dodgerblue"
        elif "watc" in f:
            inv_dict[f]["color"] = "red"
        elif "PRIMARY" in f:
            inv_dict[f]["color"] = "magenta"
            inv_dict[f]["marker"] = "o"
            inv_dict[f]["size"] = 10
            inv_dict[f]["text"] = False
        elif "AUXILIARY" in f:
            inv_dict[f]["color"] = "powderblue"
            inv_dict[f]["marker"] = "o"
            inv_dict[f]["text"] = False
            inv_dict[f]["size"] = 10

    fig, ax = setup(proj_str="Robinson")
    for key, vals in inv_dict.items():
        inv = vals.pop("inv")
        plot_inventory(inv=inv, ax=ax, **vals)
    plt.savefig("figures/ims_station_list.png", dpi=400)
    # plt.show()
