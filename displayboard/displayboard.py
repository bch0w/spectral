"""
A simple earthquake fetcher and mapper for quick glances at the most recent
earthquakes for the past X amount of days.
"""
from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import beach
from obspy import read_events, UTCDateTime

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def get_gcmt_quick(number_of_days=30, magnitude=4):
    """

    """
    cat = read_events(
        "http://www.ldeo.columbia.edu/~gcmt/projects/CMT/"
        "catalog/NEW_QUICK/qcmt.ndk"
    )
    previously_on = UTCDateTime() - (number_of_days * 60 * 60 * 24)
    cat_filt = cat.filter("time > {}".format(str(previously_on - 60)),
                          "magnitude >= {}".format(magnitude)
                          )
    return cat_filt


def event_beachball(m, event):
    """
    Plot event beachball for a given geonet moment tensor list,
    read in the from the GeoNet moment tensor csv file.

    :type m: Basemap
    :param m: basemap object
    :type event: obspy.core.event.Event
    :param event: event object which should contain focal mechanism
    """
    eventx, eventy = m(event.preferred_origin().longitude,
                       event.preferred_origin().latitude
                       )
    focal_mechanism = [
        event.preferred_focal_mechanism().moment_tensor.tensor['m_rr'],
        event.preferred_focal_mechanism().moment_tensor.tensor['m_tt'],
        event.preferred_focal_mechanism().moment_tensor.tensor['m_pp'],
        event.preferred_focal_mechanism().moment_tensor.tensor['m_rt'],
        event.preferred_focal_mechanism().moment_tensor.tensor['m_rp'],
        event.preferred_focal_mechanism().moment_tensor.tensor['m_tp']
                   ]
    b = beach(focal_mechanism, xy=(eventx, eventy), width=10, linewidth=1,
              facecolor='r')
    b.set_zorder(1000)
    ax = plt.gca()
    ax.add_collection(b)


def initiate_basemap():
    """
    set up the basemap object in the same way each time

    :type map_corners: list of floats
    :param map_corners: [lat_bot,lat_top,lon_left,lon_right]
    :rtype m: Basemap
    :return m: basemap object
    """
    # initiate map
    m = Basemap(projection='cyl', lon_0=180, lat_0=0)
    # m.drawmeridians(np.arange(0, 360, 30))
    # m.drawparallels(np.arange(-90, 90, 30))
    m.etopo()

    return m


def create_depth_colorbar(cat):
    """
    make a colorbar based on depths of events
    :param cat:
    :return:
    """
    depths = []
    for event in cat:
        depths.append(event.preferred_origin().depth * 1E-3)

    norm = mpl.colors.Normalize(vmin=min(depths), vmax=max(depths))
    cmap = cm.plasma
    colormap = cm.ScalarMappable(norm=norm, cmap=cmap)

    return colormap


def annotate_event_information(xy, event):
    """
    annotate onto plot
    :param m:
    :param event:
    :return:
    """
    template = "{origintime}\n{latitude}, {longitude}\n{magnitude}"
    plt.annotate(template.format(
        origintime=event.preferred_origin().time,
        latitude=event.preferred_origin().latitude,
        longitude=event.preferred_origin().longitude,
        magnitude=event.preferred_magnitude().mag),
        xy=xy, bbox=dict(facecolor='w', edgecolor='k'), multialignment='center',
        fontsize=6
                 )


def plot_events():
    """
    plot events on map
    :return:
    """
    m = initiate_basemap()
    # cat = get_gcmt_quick(number_of_days=30, magnitude=4)
    cat = read_events('/Users/chowbr/Documents/subduction/data/QUAKEML/displayboard_testcat.xml')
    colormap = create_depth_colorbar(cat)
    for event in cat:
        xy = m(event.preferred_origin().longitude,
               event.preferred_origin().latitude)
        if event.preferred_magnitude().mag > 6.0:
            event_beachball(m, event, colormap)
        else:
            m.scatter(xy[0], xy[1], marker='o', edgecolor='k', s=50,
                      color=colormap.to_rgba(
                          event.preferred_origin().depth*1E-3))
            annotate_event_information(xy, event)
    plt.show()


if __name__ == "__main__":
    plot_events()

