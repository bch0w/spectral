"""
Generate a publication quality map using Basemap and Matplotlib.
Capabilities for plotting stations, events, fault lines, contours,
landmarks etc. Designed to plot the North Island of New Zealand but should be
fairly general.

Requires configs to define flags and keyword arguments, located in 'mapcfgs'
To change the config file, change the import statement in main()

Tickmarks:
    https://stackoverflow.com/questions/18363987/basemap-how-to-remove-actual-lat-lon-lines-while-keeping-the-ticks-on-the-axis/62087643#62087643
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from obspy import read_events
from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import beach


class FixPointNormalize(mpl.colors.Normalize):
    """
    Taken from:
    https://stackoverflow.com/questions/40895021/python-equivalent-for-matlabs-demcmap-elevation-appropriate-colormap
    ===============================
    Inspired by
    https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib
    Subclassing Normalize to obtain a colormap with a fixpoint
    somewhere in the middle of the colormap.

    This may be useful for a `terrain` map, to set the "sea level"
    to a color in the blue/turquise range.

    col_val sets whichever value in the colormap to be the sealevel [0, 1]
    Stackoverflow decided 0.22 was a good choice

    Example Call:
    norm4 = FixPointNormalize(sealevel=0, vmax=3400)
    im4 = ax[1,0].imshow(data+1000, norm=norm4, cmap=cut_terrain_map)
    fig.colorbar(im4, ax=ax[1,0])
    """
    def __init__(self, vmin=None, vmax=None, sealevel=0, col_val=0.21875,
                 clip=False):
        # sealevel is the fix point of the colormap (in data units)
        self.sealevel = sealevel

        # col_val is the color value in the range [0,1] that should represent
        # the sealevel.
        self.col_val = col_val
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.sealevel, self.vmax], [0, self.col_val, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def myround(x, base=5, choice='near'):
    """
    Round value x to nearest base, round 'up','down' or to 'near'est base

    :type x: float
    :param x: value to be rounded
    :type base: int
    :param base: nearest integer to be rounded to
    :type choice: str
    :param choice: method of rounding, 'up', 'down' or 'near'
    :rtype roundout: int
    :return: rounded value
    """
    if choice == 'near':
        roundout = int(base * round(float(x)/base))
    elif choice == 'down':
        roundout = int(base * np.floor(float(x)/base))
    elif choice == 'up':
        roundout = int(base * np.ceil(float(x)/base))

    return roundout


def normalize_a_to_b(array, a=0, b=1):
    """
    normalize an array from a to b for e.g. plotting, maths

    :type array: list
    :param array: values to be normalized
    :type a: int
    :param a: lower bound of normalization
    :type b: int
    :param b: upper bound of normalization
    :rtype z: numpy.array
    :return z: normalized array
    """
    array = np.array(array)
    z = ((b-a) * (array-array.min()) / (array.max()-array.min())) + a

    return z


def cut_terrain_cmap(land_segments=200, sea_segments=56):
    """
    Combine the lower and upper range of the terrain colormap with a gap in the
    middle to let the coastline appear more prominently.
    inspired by
    https://stackoverflow.com/questions/31051488/combining-two-matplotlib-colormaps
    """
    colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 56))
    colors_land = plt.cm.terrain(np.linspace(0.25, 1, 200))

    # combine them and build a new colormap
    colors = np.vstack((colors_undersea, colors_land))
    cmap_segments = land_segments + sea_segments
    return mpl.colors.LinearSegmentedColormap.from_list('cuterraine', colors,
                                                        cmap_segments)

def discretize_cmap(cmap, N):
    """
    Set the number of values to some discrete number of values N
    :return:
    """
    pass
    # color_list=


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
    zorder = kwargs.get("zorder", 5)
    coastline_linewidth = kwargs.get("coastline_linewidth", 2.0)
    fill_color = kwargs.get("fill_color", "w")
    fontsize = kwargs.get("fontsize", 13)
    area_thresh = kwargs.get("area_thresh", None)
    resolution = kwargs.get("resolution", "h")
    degrees = kwargs.get("degrees", 1)

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

    m.drawrivers(linewidth=0.)
    m.drawcoastlines(linewidth=coastline_linewidth, zorder=zorder)
    if continent_color:
        m.fillcontinents(color=continent_color, lake_color=lake_color,
                         zorder=zorder)
    if fill_color:
        m.drawmapboundary(fill_color=fill_color)

    m.drawparallels(np.arange(int(map_corners["lat_min"]),
                              int(map_corners["lat_max"]), degrees),
                    labels=[1, 0, 0, 0], zorder=zorder+10,
                    dashes=[6, 900], linewidth=axis_linewidth,
                    fontsize=fontsize,
                    rotation=45)
    m.drawmeridians(np.arange(int(map_corners["lon_min"]),
                              int(map_corners["lon_max"]) + 1, degrees),
                    labels=[0, 0, 0, 1], fontsize=fontsize, zorder=zorder+10,
                    dashes=[6, 900], linewidth=axis_linewidth,
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
    scalebar_length = kwargs.get("scalebar_length", 100)

    mc = map_corners
    if loc == "upper-right":
        latscale = mc["lat_min"] + (mc["lat_max"] - mc["lat_min"]) * 0.94
        lonscale = mc["lon_min"] + (mc["lon_max"] - mc["lon_min"]) * 0.875
    if loc == "lower-right":
        # ORIGINAL
        latscale = mc["lat_min"] + (mc["lat_max"] - mc["lat_min"]) * 0.04
        lonscale = mc["lon_min"] + (mc["lon_max"] - mc["lon_min"]) * 0.9
        
        # BEACON
        # latscale = mc["lat_min"] + (mc["lat_max"] - mc["lat_min"]) * 0.1
        # lonscale = mc["lon_min"] + (mc["lon_max"] - mc["lon_min"]) * 0.8

    m.drawmapscale(lonscale, latscale, lonscale, latscale, scalebar_length,
                   yoffset=0.01 * (m.ymax-m.ymin), zorder=5000,
                   linewidth=linewidth, fontsize=fontsize
                   )


def plot_stations(m, fid, annotate=False, 
                  sta_ignore=[], net_ignore=[], **kwargs):
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

    networks = np.unique(stations[:, 1])
    # networks = {net: f"C{i}" for i, net in enumerate(networks)}

    xlist, ylist, clist = [], [], []
    for station in stations:
        sta, net, lat, lon = station
        if net in net_ignore:
            continue
        if sta in sta_ignore:
            continue
        x, y = m(float(lon), float(lat))

        xlist.append(x)
        ylist.append(y)
        if annotate and net in ["NZ"]:
            sta_ = station[0]
            ha = "center"
            if sta_ == "RD17":
                sta_ = "RD17/RD21"
                ha = "left"
            elif sta_ == "RD11":
                sta_ = "RD11/RD20"
                ha = "left"
            plt.text(x=x, y=y+3000, s=sta_, 
                     zorder=zorder+5, va="bottom", ha=ha, fontsize=12)
        if net == "NZ":
            c = "g"
            m_ = "d"
        elif net == "XX":
            c = "orange"
            m_ = "^"
        elif net == "X2":
            c = "r"
            m_ = "<"
        elif net == "ZX":
            c = "forestgreen"
            m_ = ">"
        elif net == "Z8":
            c = "yellow"
            m_ = "d"
        elif net == "SL":
            c = "g"
            m_ = "d"
        elif net == "CC":
            c = "purple"
            m_ = ">"
        else:
            c = "w"
            m_ = "v"

        clist.append(c)
        m.scatter(x, y, marker=m_, s=markersize, edgecolors="k",
                  c=c, linestyle="-", linewidth=2., zorder=zorder)

        # if color == "by_network":
        #     c.append(networks[net])
        # else:
        #     c.append(color)

    # m.scatter(xlist, ylist, marker=marker, s=markersize, edgecolors="k",
    #           c=clist, linestyle="-", linewidth=2., zorder=zorder)


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
    cbar_ticknum = kwargs.get("cbar_ticknum", 11)
    cbar_labelpad = kwargs.get("cbar_labelpad", 15)
    cbar_linewidth = kwargs.get("cbar_linewidth", 2.)
    cmap_name = kwargs.get("cmap_name", "jet_r")
    zorder = kwargs.get("zorder", 10)

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
        try:
            event_beachball(m, event, fm_type="focal_mechanism",
                            src_markercolor=cmap(normalize(depths[i])),
                            src_width=mags_norm[i], 
                            zorder=zorder + 1/mags_norm[i]
                            )
        except Exception as e:
            print(f"event {i} failed")
            continue

    # Create a colormap for the colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=normalize)
    sm.set_array([])
    cbar = plt.colorbar(sm, extend="max", shrink=cbar_shrink, aspect=10)
    # To change the number of labels
    # cbar.locator = mpl.ticker.MaxNLocator(nbins=cbar_ticknum)
    # cbar.update_ticks()
    cbar.set_label("depth (km)", rotation=270, labelpad=cbar_labelpad,
                   fontsize=cbar_fontsize)
    cbar.ax.invert_yaxis()
    cbar.ax.tick_params(axis="y", direction="in", length=0, width=3., 
                        labelsize=cbar_tickfontsize)
    cbar.outline.set_linewidth(cbar_linewidth)

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


def plot_raypaths(m, events_fid, stations_fid, pairs, **kwargs):
    """
    Plot lines connecting sources and receivers

    :type pairs: dict
    :param pairs: src rcv pairs {event_id: station_name}
    """
    event_color = kwargs.get("event_color", "r")

    cat = read_events(events_fid)
    cat_indices = [_.resource_id.id.split('/')[1] for _ in cat]

    stations = np.loadtxt(stations_fid, usecols=(0,1,2,3), dtype=str)
    sta_indices = [f"{_[1]}.{_[0]}" for _ in stations]

    plotted_stations = []

    for i, (event_id, station_name) in enumerate(pairs):
        # Plot focal mechanism
        event = cat[cat_indices.index(event_id)]
        ex, ey = m(event.preferred_origin().longitude,
                   event.preferred_origin().latitude
                   )

        # original src_width=1E4
        event_beachball(m, event, fm_type="focal_mechanism",
                        src_markercolor="r", src_width=3E4, zorder=100
                        )
        # Plot station
        station = stations[sta_indices.index(station_name)]
        _, _, lat, lon = station
        sx, sy = m(float(lon), float(lat))
        if station_name.split(".")[0] == "NZ":
            c = "c"
        else:
            c = "darkorange"

        if station_name not in plotted_stations:
            # Original s=50
            plt.scatter(sx, sy, marker="v", s=150, edgecolor="k", linewidth=2.,
                        zorder=100, color="C4")#c)
            plotted_stations.append(station_name)

        # Connect with a line (thick and colored)
        plt.plot([sx, ex], [sy, ey], color=f"C4", linestyle="-", 
                 linewidth=3., zorder=90)

        # Connect with a line (thin and black)
        # plt.plot([sx, ex], [sy, ey], color=f"k", linestyle="-", 
        #         linewidth=.75, zorder=90, alpha=0.2)

        # Find midpoint in line annotate
        if False:
            xvals = [sx, ex]
            yvals = [sy, ey]
            midpoint = ((max(xvals) - min(xvals)) / 3 + min(xvals), 
                        (max(yvals) - min(yvals)) / 3 + min(yvals))
            plt.text(x=midpoint[0], y=midpoint[1], 
                     s=f"{event_id}_{station_name}", fontsize=15, zorder=200)


        print(event_id, station_name)


def draw_domain_bounds(m, bounds, **kwargs):
    """
    Draw the domain boundaries if theyre smaller than the map boundaries
    """
    lats = [bounds["lat_min"], bounds["lat_max"], 
            bounds["lat_max"], bounds["lat_min"]]
    lons = [bounds["lon_min"], bounds["lon_min"], 
            bounds["lon_max"], bounds["lon_max"]]
    x, y = m(lons, lats)
    xy = zip(x,y)
    poly = mpl.patches.Polygon(list(xy), **kwargs)

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
    edgecolor = kwargs.get("edgecolor", "k")
    fontsize = kwargs.get("fontsize", 12)
    linewidth = kwargs.get("linewidth", 2)
    zorder = kwargs.get("zorder", 100)
    alpha = kwargs.get("alpha", 0.85)
    textborder = kwargs.get("textborder", "k")
    fontcolor = kwargs.get("fontcolor", "k")
    marker = kwargs.get("marker", "o")
    fontweight = kwargs.get("fontweight", "normal")
    annotate = kwargs.get("annotate", True)
    mark = kwargs.get("mark", True)

    for name, tup in locations.items():
        try:
            lat, lon = tup
            kwargs = {}
        except ValueError:
            lat, lon, kwargs = tup
        xy = m(lon, lat)
        if mark:
            m.scatter(xy[0], xy[1], marker=marker, s=markersize,
                      edgecolor=edgecolor, c=color, linewidth=linewidth,
                      zorder=zorder)
        if annotate:
            plt.text(xy[0] + 0.01 * (m.xmax - m.xmin),
                     xy[1] + 0.01 * (m.ymax - m.ymin),
                     s=name, fontsize=fontsize, zorder=zorder,
                     fontweight=fontweight, c=fontcolor,
                     path_effects=[pe.withStroke(linewidth=linewidth,
                                                 foreground=textborder)],
                     **kwargs
                     # bbox=dict(facecolor="w", fill=True, edgecolor="k",
                     #           linewidth=1.5, alpha=alpha, boxstyle="round")
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
    clabels = plt.clabel(cs, levels[1:], fontsize=fontsize, fmt=format_,)

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


def plot_topography(m, fid, **kwargs):
    """
    Plot topography based on an xyz file
    :param m:
    :param fid:
    :param kwargs:
    :return:
    """
    zorder = kwargs.get("zorder", 20)
    cbar_shrink = kwargs.get("cbar_shrink", 0.35)
    cbar_fontsize = kwargs.get("cbar_fontsize", 15)
    cbar_tickfontsize = kwargs.get("cbar_tickfontsize", 15)
    cbar_labelpad = kwargs.get("cbar_labelpad", 15)
    cbar_linewidth = kwargs.get("cbar_linewidth", 2.)
    land_segments = kwargs.get("land_segments", 200)
    sea_segments = kwargs.get("sea_segments", 56)
    zero_col_val = kwargs.get("zero_col_val", .21875)
    markersize = kwargs.get("markersize", 0.05)

    # Load in the topography/bathymetry file
    lon, lat, z = np.loadtxt(fid).T
    z *= 1E-3
    x, y = m(lon, lat)
    
    # Generate a segmented colormap that creates a hard boundary at coastline
    # colors (between browns and blues)
    cmap = cut_terrain_cmap(land_segments=land_segments, 
                            sea_segments=sea_segments)

    # Set 0 value to sea-level, so turquoise colors
    normalize = FixPointNormalize(sealevel=0, vmax=max(z), vmin=min(z),
                                  col_val=zero_col_val)

    # Plot that ish, small pointsize to get fine detail
    sc = m.scatter(x, y, c=z, cmap=cmap, zorder=zorder, norm=normalize, 
                   s=markersize)

    # Create a corresponding colorbar
    cbar = plt.colorbar(sc, shrink=cbar_shrink, aspect=10)  # extend="both"
    cbar.set_label("elevation [km]", rotation=270, labelpad=cbar_labelpad,
                   fontsize=cbar_fontsize)
    cbar.ax.tick_params(axis="y", direction="in", length=0, width=3.,
                        labelsize=cbar_tickfontsize)
    cbar.outline.set_linewidth(cbar_linewidth)


if __name__ == "__main__":
    from configs.raypaths_small import *

    f, m = initiate_basemap(**MAP_KWARGS)

    if BOUNDS:
        draw_domain_bounds(m, DOMAIN_BOUNDS, **DMN_KWARGS)

    if STATIONS:
        plot_stations(m, FIDS["STATIONS"], **STA_KWARGS)

    if EVENTS:
        plot_beachballs(m, FIDS["EVENTS"], **EVT_KWARGS)

    if RAYPATHS:
        plot_raypaths(m, FIDS["EVENTS"], FIDS["STATIONS"], PAIRS)

    if CITIES:
        plot_landmarks(m, CITIES_DICT, **CTY_KWARGS)

    if LANDMARKS:
        plot_landmarks(m, LANDMARKS_DICT, **LMK_KWARGS)

    if INTERFACE:
        plot_interface_contour(m, FIDS["INTERFACE"], **INT_KWARGS)

    if FAULTS:
        plot_active_faults(m, FIDS["FAULTS"], **FLT_KWARGS)

    if TOPO:
        plot_topography(m, FIDS["TOPO"], **TPO_KWARGS)

    plt.savefig(FIDS["OUTPUT"], transparent=True, dpi=DPI)
    plt.show()
