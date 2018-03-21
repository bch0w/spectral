"""module file for plotting functions
"""
import getdata

def pretty_grids(input_ax):
    """make dem grids pretty
    """
    input_ax.set_axisbelow(True)
    input_ax.tick_params(which='both',direction='in',top=True,right=True)
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
