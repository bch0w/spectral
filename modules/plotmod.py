"""module file for plotting functions
"""
import getdata
import matplotlib.pyplot as plt

def build_color_dictionary(num_of_colors,map='nipy_spectral'):
    """to keep static colors for each station for future plots. assumes 21
    stations in the array
    """
    import matplotlib.cm as mplcm
    import matplotlib.colors as colors

    num_of_colors += 1
    cm = plt.get_cmap(map)
    norm_col = colors.Normalize(vmin=0,vmax=num_of_colors-1)
    scalarMap = mplcm.ScalarMappable(norm=norm_col,cmap=cm)
    colorrange = [scalarMap.to_rgba(i) for i in range(num_of_colors)]
        
    color_dictionary = {}
    for i in range(0,num_of_colors):
        color_dictionary[i] = colorrange[i]

    return color_dictionary

def pretty_grids(input_ax):
    """make dem grids pretty
    """
    input_ax.set_axisbelow(True)
    input_ax.tick_params(which='both',
                         direction='in',
                         top=True,
                         right=True)
    input_ax.minorticks_on()
    input_ax.grid(which='minor',
                    linestyle=':',
                    linewidth='0.5',
                    color='k',
                    alpha=0.25)
    input_ax.grid(which='major',
                    linestyle='-',
                    linewidth='0.5',
                    color='k',
                    alpha=0.15)
    input_ax.ticklabel_format(style='sci',
                            axis='y',
                            scilimits=(0,0))

def align_yaxis(ax1,v1,ax2,v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1
    """
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)


def event_beachball(eventid,fig):
    """plt event beachball on figure object
    """
    MT = getdata.get_moment_tensor(eventid)
    eventx,eventy = fig.bmap(MT['Longitude'],MT['Latitude'])
    FM = [MT['strike1'],MT['dip1'],MT['rake1']]

    b = beach(FM,xy=(eventx,eventy),width=3E4,linewidth=1,facecolor='r')
    b.set_zorder(10)
    ax = plt.gca()
    ax.add_collection(b)
    plt.annotate("{}".format(eventid),
                    xy=(eventx,eventy),
                    xytext=(eventx,eventy),
                    fontsize=7,
                    zorder=200,
                    weight='bold')


def compare_beachballs(event_id):
    """plot beachballs to compare GCMT and converted GEONET MT
    """
    from obspy.imaging.beachball import beachball
    from getdata import get_moment_tensor
    import matplotlib.pyplot as plt

    MT = get_moment_tensor(event_id=event_id)
    geonet_xyz = [MT['Mxx'],MT['Myy'],MT['Mzz'],MT['Mxy'],MT['Mxz'],MT['Myz']]
    geonet_rtp = mt_transform(geonet_xyz,method='xyz2rtp')

    gcmt = getdata.get_GCMT_solution(event_id=event_id)
    GCMT = mt_from_event(gcmt)
    gcmt_rtp = [GCMT['Mrr'],GCMT['Mtt'],GCMT['Mpp'],
                GCMT['Mrt'],GCMT['Mrp'],GCMT['Mtp']]
    gcmt_xyz = mt_transform(geonet_rtp,method='rtp2xyz')

    bp = pathnames()['kupeplots'] + 'beachballs/{e}_{o}_{c}.png'
    import ipdb;ipdb.set_trace()
    geonet_xyz_b = beachball(geonet_xyz,outfile=bp.format(e=event_id,
                                                          o='geonet',
                                                          c='xyz',
                                                          ))
    geonet_rtp_b = beachball(geonet_rtp,outfile=bp.format(e=event_id,
                                                          o='geonet',
                                                          c='rtp',
                                                          ))
    gcmt_rtp_b = beachball(gcmt_rtp,outfile=bp.format(e=event_id,
                                                      o='gcmt',
                                                      c='rtp',
                                                      ))
    gcmt_xyz_b = beachball(gcmt_xyz,outfile=bp.format(e=event_id,
                                                      o='gcmt',
                                                      c='xyz',
                                                      ))
    plt.close("all")


def plot_event_station(inv,cat):
    """plot event and station along with line connection and station-receiver
    distance for quick visualization
    """
    from obspy.geodetics.base import gps2dist_azimuth
    # plot station on a map
    fig = inv.plot(marker='v',
            projection="local",
            resolution = "i",
            size=50,
            show=False,
            continent_fill_color="white",
            water_fill_color="white",
            color_per_network=True,
            label=False)

    # convert station and event information to map coords
    station_lat = inv[0][0].latitude
    station_lon = inv[0][0].longitude
    station_x,station_y = fig.bmap(station_lon,station_lat)
    event_lat = cat[0].origins[0].latitude
    event_lon = cat[0].origins[0].longitude
    event_x,event_y = fig.bmap(event_lon,event_lat)

    # annotate station code
    # !!!!!!!!!!!!!!!!!!!!!!!!! only for waveform_by_event.py
    # for INV in inv[0]:
    #     station_lat = INV.latitude
    #     station_lon = INV.longitude
    #     station_x,station_y = fig.bmap(station_lon,station_lat)
    #     stationcode = INV.code
    #     dist_azi = gps2dist_azimuth(event_lat,event_lon,station_lat,station_lon)
    #     epi_dist = round(dist_azi[0]*1E-3,2)
    #     BAz = round(dist_azi[2],2)
    #     anno_template = "{}\n{}km".format(stationcode,epi_dist)
    #     plt.annotate(anno_template,
    #                 xy=(station_x,station_y),
    #                 xytext=(station_x,station_y),
    #                 fontsize=10,
    #                 weight='bold',
    #                 zorder=100)
    # !!!!!!!!!!!!!!!!!!!!!!!!!

    stationcode = inv[0][0].code
    plt.annotate(stationcode,
                xy=(station_x,station_y),
                xytext=(station_x,station_y),
                fontsize=10,
                weight='bold',
                zorder=100)

    # plot event
    scatter = fig.bmap.scatter(event_x,event_y,
                                marker='o',
                                s=150,
                                zorder=50,
                                edgecolor='k')
    magnitude = round(cat[0].magnitudes[0].mag,2)
    plt.annotate(magnitude,
                xy=(event_x,event_y),
                xytext=(event_x,event_y),
                fontsize=10,
                weight='bold',
                zorder=100)

    # connect station and event by arrow
    dist_azi = gps2dist_azimuth(event_lat,event_lon,station_lat,station_lon)
    epi_dist = round(dist_azi[0]*1E-3,2)
    BAz = round(dist_azi[2],2)
    plt.title("Epicentral Distance: {} km | BAz: {}".format(epi_dist,BAz))


    plt.show()
