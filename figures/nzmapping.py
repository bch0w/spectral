"""
Generate a publication quality map using Basemap and Matplotlib.
Capabilities for plotting stations, events, fault lines, contours,
landmarks etc. Designed to plot the North Island of New Zealand but should be
fairly general.

Requires configs to define flags and keyword arguments, located in 'mapcfgs'
To change the config file, change the import statement in main()
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from obspy import read_events
from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import beach
from pyatoa.utils.calculate import normalize_a_to_b, myround


def initiate_basemap(map_corners, scalebar=True, **kwargs):
    """
    set up the basemap object in the same way each time

    :type map_corners: dict of floats
    :param map_corners: {lat_min,lat_max,lon_min,lon_max}
    :type scalebar: bool
    :param scalebar: add a scalebar to the map
    :rtype m: Basemap
    :return m: basemap object
    """
    figsize = kwargs.get("figsize", (8, 10))
    dpi = kwargs.get("dpi", 100)
    axis_linewidth = kwargs.get("axis_linewidth", 2)
    continent_color = kwargs.get("continent_color", "w")
    lake_color = kwargs.get("lake_color", "w")
    coastline_zorder = kwargs.get("coastline_zorder", 5)
    coastline_linewidth = kwargs.get("coastline_linewidth", 2.0)
    fill_color = kwargs.get("fill_color", "w")
    fontsize = kwargs.get("fontsize", 13)
    area_thresh = kwargs.get("area_thresh", None)
    resolution = kwargs.get("resolution", "h")

    f = plt.figure(figsize=figsize, dpi=dpi)
    ax = f.add_subplot(111)

    # Initiate map and draw in style
    m = Basemap(projection="stere", resolution=resolution, rsphere=6371200,
                area_thresh=area_thresh,
                lat_0=(map_corners["lat_min"] + map_corners["lat_max"]) / 2,
                lon_0=(map_corners["lon_min"] + map_corners["lon_max"]) / 2,
                llcrnrlat=map_corners["lat_min"],
                urcrnrlat=map_corners["lat_max"],
                llcrnrlon=map_corners["lon_min"],
                urcrnrlon=map_corners["lon_max"],
                )

    m.drawcoastlines(linewidth=coastline_linewidth, zorder=coastline_zorder)
    m.fillcontinents(color=continent_color, lake_color=lake_color)
    m.drawrivers(linewidth=0.)
    m.drawmapboundary(fill_color=fill_color)
    m.drawparallels(np.arange(int(map_corners["lat_min"]),
                              int(map_corners["lat_max"]), 1),
                    labels=[1, 0, 0, 0], linewidth=0, fontsize=fontsize,
                    rotation=45)
    m.drawmeridians(np.arange(int(map_corners["lon_min"]),
                              int(map_corners["lon_max"]) + 1, 1),
                    labels=[0, 0, 0, 1], linewidth=0, fontsize=fontsize,
                    rotation=45)

    if scalebar:
        place_scalebar(m, map_corners, **kwargs)

    ax = plt.gca()
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(axis_linewidth)

    f.tight_layout()

    return f, m


def place_scalebar(m, map_corners, **kwargs):
    """
    Put the scale bar in a corner at a reasonable distance from each edge

    Handy reminder for moving the scalebar around:
        latitude is up, down
        longitude is right, left

    :type m: Basemap
    :param m: basemap object
    :type map_corners: dict of floats
    :param map_corners: [lat_bot,lat_top,lon_left,lon_right]
    :type loc: str
    :param loc: location of scalebar, "upper-right" or "lower-right"
    """
    loc = kwargs.get("scalebar_location", "upper-right")
    fontsize = kwargs.get("scalebar_fontsize", 13)
    linewidth = kwargs.get("scalebar_linewidth", 2)

    mc = map_corners
    if loc == "upper-right":
        latscale = mc["lat_min"] + (mc["lat_max"] - mc["lat_min"]) * 0.94
        lonscale = mc["lon_min"] + (mc["lon_max"] - mc["lon_min"]) * 0.875
    if loc == "lower-right":
        latscale = mc["lat_min"] + (mc["lat_max"] - mc["lat_min"]) * 0.04
        lonscale = mc["lon_min"] + (mc["lon_max"] - mc["lon_min"]) * 0.9

    m.drawmapscale(lonscale, latscale, lonscale, latscale, 100,
                   yoffset=0.01 * (m.ymax-m.ymin), zorder=5000,
                   linewidth=linewidth, fontsize=fontsize
                   )


def plot_stations(m, fid, sta_ignore=[], net_ignore=[], **kwargs):
    """
    Get station info from a Specfem STATION file
    
    :type m: basemap
    :param m: map for converting coordinates
    :type station_info: np.array
    :param station_info: array from reading in the station file using numpy
    """
    markersize = kwargs.get("markersize", 100)
    color = kwargs.get("color", "w")
    fontsize = kwargs.get("fontsize", 10)
    zorder = kwargs.get("zorder", 100)
    edgecolor = kwargs.get("station_edgecolor", "k")
    marker = kwargs.get("station_marker", "v")

    stations = np.loadtxt(fid, usecols=(0,1,2,3), dtype=str) 

    x, y = [], []
    for station in stations:
        sta, net, lat, lon = station
        if net in net_ignore:
            continue
        if sta in sta_ignore:
            continue
        xy = m(float(lon), float(lat))

        x.append(xy[0])
        y.append(xy[1])

    m.scatter(x, y, marker=marker, s=markersize, edgecolors=edgecolor,
              c=color, linestyle="-", linewidth=2., zorder=zorder)


def event_beachball(m, event, fm_type="focal_mechanism", zorder=90, **kwargs):
    """
    Plot event beachball for a given moment tensor attribute from event object.
    Note:
    if strike_dip_rake chosen, nodal plane 1 is used for focal mechanism,
    assuming that is the preferred plane
    :type m: Basemap
    :param m: basemap object
    :type fm_type: str
    :param fm_type: focal mech. type, "focal_mechanism" or "strike_dip_rake"
    :type event: obspy.core.event.Event
    :param event: event object which should contain focal mechanism
    """
    width = kwargs.get("src_width", 2.6E4)

    src_marker = kwargs.get("src_marker", "o")
    src_markercolor = kwargs.get("src_markercolor", "r")
    src_markersize = kwargs.get("src_markersize", 105)
    src_linewidth = kwargs.get("src_linewidth", 1.75)

    eventx, eventy = m(event.preferred_origin().longitude,
                       event.preferred_origin().latitude
                       )

    # No focal mechanism? Just plot a ploint, same as connect_source_receiver()
    if not hasattr(event, "focal_mechanisms"):
        m.scatter(eventx, eventy, marker=src_marker, color=src_markercolor,
                  edgecolor="k", s=src_markersize, linewidth=src_linewidth,
                  zorder=90)

    else:
        if fm_type == "focal_mechanism":
            fm = event.focal_mechanisms[0].moment_tensor.tensor or \
                 event.preferred_focal_mechanism().moment_tensor.tensor
            beach_input = [fm["m_rr"], fm["m_tt"], fm["m_pp"],
                           fm["m_rt"], fm["m_rp"], fm["m_tp"]
                           ]
        elif fm_type == "strike_dip_rake":
            nod_plane = event.focal_mechanisms[0].nodal_planes or \
                        event.preferred_focal_mechanism().nodal_planes
            # try determine the preferred nodal plane, default to 1
            try:
                sdr = nod_plane[f"nodal_plane_{nod_plane.preferred_plane}"]
            except AttributeError:
                sdr = nod_plane.nodal_plane_1
            beach_input = [sdr.strike, sdr.dip, sdr.rake]

        b = beach(beach_input, xy=(eventx, eventy), width=width,
                  linewidth=src_linewidth, facecolor=src_markercolor
                  )
        b.set_zorder(zorder)
        ax = plt.gca()
        ax.add_collection(b)


def plot_beachballs(m, fid, norm_a, norm_b, mag_scale=None,
                    mag_legend=True,**kwargs):
    """
    Plot moment tensors as beachballs, color by depth if requested

    :type m: basemap
    :param m: map to plot beacbhalls on
    :type cat: obspy.Catalog
    :param cat: catalog with events
    :type depths: list
    :param depths: depths in km of events
    :type mags: list
    :param mags: magnitudes of events
    :type cmap: matplotlib.colors.ListedColorMap
    :param cmap: chosen colormap for event depths
    :type normalize: matplotlib.colors.BoundaryNorm
    :param normalize: normalized discrete colormap for colorbar
    """
    cbar_shrink = kwargs.get("cbar_shrink", 0.35)
    cbar_fontsize = kwargs.get("cbar_fontsize", 15)
    cbar_tickfontsize = kwargs.get("cbar_tickfontsize", 15)
    cbar_labelpad = kwargs.get("cbar_labelpad", 15)
    cmap_name = kwargs.get("cmap_name", "jet_r")

    cmap = getattr(plt.cm, cmap_name)

    cat = read_events(fid)

    # First get some array info for scale bars and colormap
    print(f"{len(cat)} events in catalog")
    ignore, depths, mags, depths = [], [], [], []
    for i, event in enumerate(cat):
        try:
            depth = event.preferred_origin().depth * 1E-3
            mag = event.preferred_magnitude().mag
        except AttributeError:
            depth = event.origins[0].depth * 1E-3
            mag = event.magnitudes[0].mag

        # if (depth > 60) or (mag < 4.4):
        #     ignore.append(i)
        #     depths.append(np.nan)
        #     mags.append(np.nan)
        #     continue

        depths.append(depth)
        mags.append(mag)

    print(f"{len(ignore)} events ignored")
   
    # Normalize the colormap, create a discrete colorbar
    max_depth = myround(max(depths), 5, "up")
    normalize = mpl.colors.BoundaryNorm(range(0, max_depth, 5), cmap.N)

    # Normalize magnitudes for sense of event size
    mags = np.array(mags)
    mags[np.isnan(mags)] = np.nanmin(mags)
    mags_norm = normalize_a_to_b(mags, a=norm_a, b=norm_b)
    
    # Plot beachballs for each event
    for i, event in enumerate(cat):
        if i in ignore:
            continue
        event_beachball(m, event, fm_type="focal_mechanism",
                        src_markercolor=cmap(normalize(depths[i])),
                        src_width=mags_norm[i], zorder=90 + 1/mags_norm[i]
                        )

    # Create a colormap for the colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=normalize)
    sm.set_array([])
    cbar = plt.colorbar(sm, extend="max", shrink=cbar_shrink, aspect=10)
    cbar.set_label("depth (km)", rotation=270, labelpad=cbar_labelpad, 
                   fontsize=cbar_fontsize)
    cbar.ax.invert_yaxis()
    cbar.ax.tick_params(axis="y", direction="in", length=0, width=3., 
                        labelsize=cbar_tickfontsize)

    # Find the bounds of the magnitude range for scale
    if mag_scale is None:
        min_mag = myround(min(mags), 0.25, "down")
        print(f"minimum magnitude: {min(mags)} -> {min_mag}")
        max_mag = myround(max(mags), 0.25, "up")
        print(f"maximum magnitude: {max(mags)} -> {max_mag}")
        mag_scale = np.arange(min_mag, max_mag + 0.25, 1)

    mag_scale_norm = normalize_a_to_b(mag_scale, a=norm_a, b=norm_b)

    # Create magnitude scale
    if mag_legend:
        dx = (m.urcrnrx - m.llcrnrx) + m.llcrnrx
        x_val = 0.9 * dx
        y_val = 0.1 * (m.urcrnry - m.llcrnry) + m.llcrnry

        for i, ms in enumerate(mag_scale_norm):
            b = beach([0, 90, 0], xy=(x_val, y_val), width=ms, facecolor="k")
            b.set_zorder(90)
            ax = plt.gca()
            ax.add_collection(b)
            mag_val = int(mag_scale[i])
            plt.text(x_val + 0.03 * dx, y_val, f"M{mag_val}", fontsize=15)

            y_val += 0.05 * (m.urcrnry - m.llcrnry) + m.llcrnry


def plot_raypaths(m, cat_fid, sta_fid, **kwargs):
    """
    Plot lines connecting sources and receivers
    """
    event_color = kwargs.get("event_color", "r")
    station_color = kwargs.get("station_color", "g")
    ray_color = kwargs.get("station_color", "k")

    # Get event lat/lon values
    cat = read_events(fid)

    # First get some array info for scale bars and colormap
    print(f"{len(cat)} events in catalog")
    eventx, eventy = [], []
    for i, event in enumerate(cat):
        try:
            origin = event.preferred_origin()
            mag = event.preferred_magnitude().mag
        except AttributeError:
            origin = event.origins[0]
            mag = event.magnitudes[0].mag

        if (origin.depth * 1E-3 > 60) or (mag < 4.4):
            continue

        x_, y_ = m(origin.longitude, origin.latitude)
        eventx.append(x)
        eventy.append(y)

    # Get station lat lon values
    stations = np.loadtxt(fid, usecols=(0,1,2,3), dtype=str)
    stax, stay = [], []
    for station in stations:
        sta, net, lat, lon = station

        x_, y_ = m(float(lon), float(lat))
        stax.append(x_)
        stay.append(y_)

    # Plot event markers, station markers and connecting line
    stations_plotted, events_plotted = [], []
    for ex, ey in zip(eventx, eventy):
        for sx, sy in zip(stax, stay):
            # Plot a marker for each event and station
            if (ex, ey) not in events_plotted:
                plt.scatter(ex, ey, marker="o", c=event_color, edgecolors="k")
                events_plotted.append((ex, ey))
            if (sx, sy) not in stations_plotted:
                plt.scatter(sx, sy, marker="v", c=station_color, edgecolors="k",
                            s=25, zorder=100)
                stations_plotted.append((sx, sy))

        # Connect source and receiver with a line
        plt.plot([sx, sy], [ex, ey], color=ray_color, linestyle="-",
                 alpha=0.1, zorder=50)


def draw_domain_bounds(m, bounds):
    """
    Draw the domain boundaries if theyre smaller than the map boundaries
    """
    lats = [bounds["lat_min"], bounds["lat_max"], 
            bounds["lat_max"], bounds["lat_min"]]
    lons = [bounds["lon_min"], bounds["lon_min"], 
            bounds["lon_max"], bounds["lon_max"]]
    x, y = m(lons, lats)
    xy = zip(x,y)
    poly = mpl.patches.Polygon(list(xy), facecolor="None", edgecolor="k", 
                               linewidth=1.5, linestyle="--", alpha=1)

    plt.gca().add_patch(poly)


def plot_landmarks(m, locations, **kwargs):
    """
    Plot major population centers, landmarks, etc, based on user input locations

    :type m: basemap
    :param m: map to plot cities on
    :type locations: dict
    :param locations: dict of tuple which should be in the format
        {"NAME": (LAT, LON, LOC_BOOL)...}
        where LOC_BOOL is either a 1 for a city, 0 for a location
    """
    markersize = kwargs.get("markersize", 80)
    color = kwargs.get("color", "w")
    fontsize = kwargs.get("fontsize", 12)
    zorder = kwargs.get("zorder", 100)
    alpha = kwargs.get("alpha", 0.85)
    marker = kwargs.get("marker", "o")

    for name, tup in locations.items():
        lat, lon, locbool = tup
        xy = m(lon, lat)
        if locbool:
            plt.text(xy[0] + 0.01 * (m.xmax - m.xmin),
                     xy[1] + 0.01 * (m.ymax - m.ymin),
                     s=name, fontsize=fontsize, zorder=zorder,
                     fontweight="normal", bbox=dict(facecolor="w", fill=True,
                                                    edgecolor="k",
                                                    linewidth=1.5, alpha=alpha,
                                                    boxstyle="square")
                     )
        else:
            m.scatter(xy[0], xy[1], marker=marker, s=markersize, edgecolor="k",
                      c=color, linewidth=2, zorder=zorder)
            plt.text(xy[0] + 0.01 * (m.xmax - m.xmin),
                     xy[1] + 0.01 * (m.ymax - m.ymin),
                     s=name, fontsize=fontsize, zorder=zorder,
                     fontweight="normal", bbox=dict(facecolor="w", fill=True,
                                                    edgecolor="k",
                                                    linewidth=1.5, alpha=alpha,
                                                    boxstyle="round")
                     )

def plot_interface_contour(m, fid, **kwargs):
    """
    Plot the plate interface as contours, inline labels for depths

    :type m: basemap
    :param m: map convert coordinates
    :type plate: np.array
    :param plate: info read in from plate interface npy file
    :type xi, yi, zi: np.array
    :param xi, yi, zi: gridded data for use in contour
    """
    zorder = kwargs.get("zorder", 100)
    levels = kwargs.get("levels", [0, 5, 10, 15, 25, 50, 100, 200, 300])
    color = kwargs.get("color", "k")
    linestyle = kwargs.get("linestyle", "dotted")
    linewidth = kwargs.get("linewidth", 3)
    alpha = kwargs.get("alpha", 0.5)
    fontsize = kwargs.get("fontsize", 12)
    format_ = kwargs.get("format", "%.0f km")

    plate = np.load(fid)

    lats = plate[:, 0]
    lons = plate[:, 1]
    depths = plate[:, 2]
    depths = [-1 * _ for _ in depths]

    x, y = m(lons, lats)
    xi = np.linspace(min(x), max(x), 1000)
    yi = np.linspace(min(y), max(y), 1000)

    xi, yi = np.meshgrid(xi, yi)
    zi = griddata((x, y), depths, (xi, yi))

    cs = m.contour(xi, yi, zi, levels, vmin=0, colors=color,
                   linestyles=linestyle, alpha=alpha, zorder=zorder,
                   linewidths=linewidth)
    clabels = plt.clabel(cs, levels[1:], fontsize=fontsize, fmt=format_)

    [txt.set_bbox(dict(facecolor="white", edgecolor="None", pad=0))
     for txt in clabels]


def plot_active_faults(m, fid, **kwargs):
    """
    Plot GNS Active Faults
    """
    faults, lat, lon = np.loadtxt(fid).T
    x, y = m(lon, lat)

    for i in np.unique(faults):
        idx = np.where(faults == i)
        m.plot(x[idx], y[idx], **kwargs)


if __name__ == "__main__":
    from mapcfgs.srcrcv import *

    f, m = initiate_basemap(**MAP_KWARGS)

    if BOUNDS:
        draw_domain_bounds(m, DOMAIN_BOUNDS)

    if STATIONS:
        plot_stations(m, FIDS["STATIONS"], **STA_KWARGS)

    if EVENTS:
        plot_beachballs(m, FIDS["EVENTS"], **EVT_KWARGS)

    if LANDMARKS:
        plot_landmarks(m, LANDMARKS, **LMK_KWARGS)

    if INTERFACE:
        plot_interface_contour(m, FIDS["INTERFACE"], **INT_KWARGS)

    if FAULTS:
        plot_active_faults(m, FIDS["FAULTS"], **FLT_KWARGS)

    plt.savefig(FIDS["OUTPUT"])
    plt.show()
