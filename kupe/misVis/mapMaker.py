"""misfit visualization tool to be called through adjointBuilder
produces a basemap with beachball and all available stations as well as the
relevant station highlighted. important information annotated (e.g. 
misift information, distance, BAz etc.)
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from obspy.imaging.beachball import beach
from obspy.geodetics import gps2dist_azimuth

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.linewidth'] = 2


def trace_trench(m):
    """trace the hikurangi trench on a basemap object 'm'
    """
    trenchcoordspath = pathnames()['data'] + \
                                        'KUPEDATA/hikurangiTrenchCoords.npz'
    trenchcoords = np.load(trenchcoordspath)
    lats = trenchcoords['LAT']
    lons = trenchcoords['LON']
    x,y = m(lons,lats)

    # interpolate points to make a smoother curve
    xprime = np.flip(x,axis=0)
    yprime = np.flip(y,axis=0)
    xprimenew = np.linspace(x.min(),x.max(),100)
    yprimenew = np.interp(xprimenew,xprime,yprime)

    m.plot(xprimenew,yprimenew,'--',linewidth=1.25,color='k',zorder=2)    

def event_beachball(m,MT):
    """plot event beachball on basemap 'm' object for a given geonet event_id
    """
    eventx,eventy = m(MT['Longitude'],MT['Latitude'])
    FM = [MT['strike2'],MT['dip2'],MT['rake2']]

    b = beach(FM,xy=(eventx,eventy),
                 width=2.5E4,
                 linewidth=1,
                 facecolor='r')
    b.set_zorder(1000)
    ax = plt.gca()
    ax.add_collection(b)


def source_receiver(event,inv):
    """determine source receiver parameters such as great circle distance,
    backazimuth etc.
    source = A, receiver = B
    """
    MT = get_moment_tensor(event_id=event_id)
    event_lat,event_lon = event.origins[0].latitude,event.origins[0].longitude
    sta_lat,sta_lon = inv[0][0].latitude,inv[0][0].longitude
    GCDist,Az,BAz = gps2dist_azimuth(event_lat,event_lon,sta_lat,sta_lon)
    
    # get magnitude M
    for mag in event.magnitudes:
        if m.magnitude_type == "M":
            magnitude = m.mag
    
    # dictionary output for use in annotations
    srcrcvdict = {"distance":GCDist,
                  "backazimuth":BAz,
                  "date":event.origins[0].time,
                  "depth":event.origins[0].depth*1E-3,
                  "magnitude":magnitude,
                  }
    
    return srcrcvdict, MT
    

def initiate_basemap(map_corners=[-50,-32.5,165,180],figsize=(10,9.4),dpi=100):
    """set up local map of NZ to be filled
    default map corners give a rough box around new zealand
    etopo_else_flat_map: bool to control continent fill with either low res topo
    or solid line representation
    """
    # flat map configuration
    continent_color = 'w'
    lake_color = 'w'

    # initiate map
    fig = plt.figure(figsize=figsize,dpi=dpi)
    m = Basemap(projection = 'stere',
                resolution = 'h',
                rsphere = 6371200,
                lat_0 = np.mean(map_corners[:2]),
                lon_0 = np.mean(map_corners[2:]),
                llcrnrlat = map_corners[0],
                llcrnrlon = map_corners[2],
                urcrnrlat = map_corners[1],
                urcrnrlon = map_corners[3],
                )

    m.fillcontinents(color=continent_color,lake_color=lake_color)
    m.drawcoastlines(linewidth=1.5)
    trace_trench(m)

    return fig, m

def populate_basemap(m,lats,lons,names=None):
    """fill map with latitude/longitude pairs, i.e. stations, events
    """
    if names is None: names = []

    X,Y = m(lons,lats)
    scatter = m.scatter(X,Y,
                        marker='v',
                        color=w,
                        edgecolor='k',
                        s=50,
                        zorder=5)

    if len(names) != 0:
        for n_,x_,y_ in zip(names,X,Y):
            plt.annotate(n_,xy=(x_,y_),xytext=(x_,y_),zorder=6,fontsize=8.5)


def generate_map(event,inv,corners=[-42.5007,-36.9488,172.9998,179.5077],
                                                                      **kwargs):
    """initiate and populate a basemap object for New Zealands north island.
    Functionality to manually ignore stations based on user quality control
    Takes station coordinates and coloring from npz files
    Choice to annotate two stations which correspond to waveforms
    Calls beachball and trench tracer to populate basemap
    :type corners: list of floats
    :param corners: values for map corners to set bounds
     e.g. [lat_bot,lat_top,lon_left,lon_right]
    """
    PD = kwargs['PD']
    
    f,m = initiate_basemap(map_corners=corners)
    srcrcvdict,MT = source_receiver(event,inv,PD)
    event_beachball(m,MT)
    
    stationfile = pathnames()['data'] + 'STATIONXML/GEONET_AND_FATHOM.npz'
    stationlist = np.load(stationfile)
    
    populate_basemap(m,stationlist['LAT'],stationlist['LON'])
    
    # final plot additions
    m.drawmapscale(179,-41.75, 179,-41.75, 100,yoffset=0.01*(m.ymax-m.ymin))
    # plt.title(PD['event_id'])
    
    
_test_generate_map():
    
        
if __name__ == "__main__":
    generate_map()