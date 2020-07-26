"""
Make a map with events and receivers
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from obspy import read_events
from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import beach
from pyatoa.visuals.map_maker import event_beachball
from pyatoa.utils.calculate import normalize_a_to_b, myround


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
    :param loc: location of scalebar, 'upper-right' or 'lower-right'
    """
    loc = kwargs.get("scalebar_location", "upper-right")
    fontsize = kwargs.get("scalebar_fontsize", 13)
    linewidth = kwargs.get("scalebar_linewidth", 2)

    mc = map_corners
    if loc == "upper-right":
        latscale = mc['lat_min'] + (mc['lat_max'] - mc['lat_min']) * 0.94
        lonscale = mc['lon_min'] + (mc['lon_max'] - mc['lon_min']) * 0.875
    if loc == "lower-right":
        latscale = mc['lat_min'] + (mc['lat_max'] - mc['lat_min']) * 0.04
        lonscale = mc['lon_min'] + (mc['lon_max'] - mc['lon_min']) * 0.9
    m.drawmapscale(lonscale, latscale, lonscale, latscale, 100,
                   yoffset=0.01 * (m.ymax-m.ymin), zorder=5000, 
                   linewidth=linewidth, fontsize=fontsize
                   )


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
    continent_color = kwargs.get("continent_color", "w")
    lake_color = kwargs.get("lake_color", "w")
    coastline_zorder = kwargs.get("coastline_zorder", 5)
    coastline_linewidth = kwargs.get("coastline_linewidth", 2.0)
    fill_color = kwargs.get("fill_color", "w")
    fontsize = kwargs.get("fontsize", 13)

    # Initiate map and draw in style
    m = Basemap(projection='stere', resolution='h', rsphere=6371200,
                area_thresh=1000,
                lat_0=(map_corners['lat_min'] + map_corners['lat_max'])/2,
                lon_0=(map_corners['lon_min'] + map_corners['lon_max'])/2,
                llcrnrlat=map_corners['lat_min'],
                urcrnrlat=map_corners['lat_max'],
                llcrnrlon=map_corners['lon_min'],
                urcrnrlon=map_corners['lon_max'],
                )
    m.drawcoastlines(linewidth=coastline_linewidth, zorder=coastline_zorder)
    m.fillcontinents(color=continent_color, lake_color=lake_color)
    m.drawmapboundary(fill_color=fill_color)
    m.drawparallels(np.arange(int(map_corners['lat_min']),
                              int(map_corners['lat_max']), 1),
                    labels=[1, 0, 0, 0], linewidth=0, fontsize=fontsize,
                    rotation=45)
    m.drawmeridians(np.arange(int(map_corners['lon_min']),
                              int(map_corners['lon_max'])+1, 1),
                    labels=[0, 0, 0, 1], linewidth=0, fontsize=fontsize,
                    rotation=45)

    if scalebar:
        place_scalebar(m, map_corners, **kwargs)
        
    ax = plt.gca()
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(3.)

    return m


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

    m.scatter(x, y, marker=marker, s=markersize, edgecolors='k', 
              c=color, linestyle='-', linewidth=2., zorder=zorder)
    
    
def plot_beachballs(m, fid, cmap, norm_a, norm_b, mag_scale=None, **kwargs):
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

        if (depth > 60) or (mag < 4.4):
            ignore.append(i)
            depths.append(np.nan)
            mags.append(np.nan)
            continue

        depths.append(depth)
        mags.append(mag)

    print(f"{len(ignore)} events ignored")
   
    # Normalize the colormap, create a discrete colorbar
    max_depth = myround(max(depths), 5, 'up')
    normalize = mpl.colors.BoundaryNorm(range(0, max_depth, 5), cmap.N)

    # Normalize magnitudes for sense of event size
    mags = np.array(mags)
    mags[np.isnan(mags)] = np.nanmin(mags)
    mags_norm = normalize_a_to_b(mags, a=norm_a, b=norm_b)
    
    # Plot beachballs for each event
    for i, event in enumerate(cat):
        if i in ignore:
            continue
        event_beachball(m, event, fm_type="strike_dip_rake",
                        src_markercolor=cmap(normalize(depths[i])),
                        src_width=mags_norm[i], zorder=90 + 1/mags_norm[i]
                        )

    # Create a colormap for the colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=normalize)
    sm.set_array([])
    cbar = plt.colorbar(sm, extend='max', shrink=cbar_shrink, aspect=10)
    cbar.set_label("depth (km)", rotation=270, labelpad=cbar_labelpad, 
                   fontsize=cbar_fontsize)
    cbar.ax.invert_yaxis()
    cbar.ax.tick_params(axis='y', direction='in', length=0, width=3., 
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
    x_val = 0.9 * (m.urcrnrx - m.llcrnrx) + m.llcrnrx
    y_val = 0.1 * (m.urcrnry - m.llcrnry) + m.llcrnry

    for ms in mag_scale_norm:
        b = beach([0, 90, 0], xy=(x_val, y_val), width=ms, facecolor="k")
        b.set_zorder(90)
        ax = plt.gca()
        ax.add_collection(b)

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
    poly = mpl.patches.Polygon(list(xy), facecolor='None', edgecolor="k", 
                               linewidth=1.5, linestyle='--', alpha=1)

    plt.gca().add_patch(poly)


if __name__ == "__main__":
    # Paths set here
    sta_fid = ("/Users/Chow/Documents/academic/vuw/data/specfem/stations/"
               "STATIONS_NZ_NORTH")
    cat_fid = ("/Users/Chow/Documents/academic/vuw/data/events/"
               "nz_north_cat_253.xml")

    # Figure kwargs
    figsize = (10, 12)
    dpi = 100

    # Map kwargs
    map_corners={'lat_min': -42.75, 'lat_max': -36.75,  
                 'lon_min': 172.5, 'lon_max': 178.8}
    domain_bounds={'lat_min': -42.5, 'lat_max': -37.0,  
                   'lon_min': 173.0, 'lon_max': 178.5}

    coastline_linewidth = 2
    continent_color= "whitesmoke"
    lake_color = "azure"
    fill_color = "azure"
    scalebar_location = "lower-right"
    scalebar_fontsize = 18.
    scalebar_linewidth = 2.75
    fontsize = 22.5

    # Station kwargs
    station_markersize = 90
    station_zorder = 100
    station_color = "w"
    station_marker = "v"
    station_edgecolor = "k"

    # Event kwargs
    beachball_cmap = plt.cm.jet_r
    beachball_norm_a = 1.25E4
    beachball_norm_b = 2.75E4
    beachball_linewidths = 1.6
    beachball_zorder = station_zorder - 2
    cbar_shrink = 0.2
    cbar_fontsize = 15
    cbar_tickfontsize = 15
    cbar_labelpad = 17.5


    # Plotting starts here
    f = plt.figure(figsize=figsize, dpi=dpi)
    ax = f.add_subplot(111)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(2)

    m = initiate_basemap(map_corners, continent_color=continent_color,
                         lake_color=lake_color, fill_color=fill_color,
                         coastline_linewidth=coastline_linewidth, scalebar=True,
                         scalebar_location=scalebar_location, fontsize=fontsize,
                         scalebar_fontsize=scalebar_fontsize, 
                         scalebar_linewidth=scalebar_linewidth
                         )
    draw_domain_bounds(m, domain_bounds)

    plot_stations(m, sta_fid, markersize=station_markersize, 
                  zorder=station_zorder, color=station_color, 
                  station_edgecolor=station_edgecolor,  
                  station_marker=station_marker
                  )

    plot_beachballs(m, cat_fid, cmap=beachball_cmap, norm_a=beachball_norm_a, 
                    norm_b=beachball_norm_b, linewidth=beachball_linewidths, 
                    zorder=beachball_zorder, cbar_shrink=cbar_shrink
                    )

    f.tight_layout()
    plt.savefig("big_map_raw.png")
    plt.show()
