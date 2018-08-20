"""misfit visualization tool to be called through adjointBuilder
produces a basemap with beachball and all available stations as well as the
relevant station highlighted. important information annotated (e.g.
misift information, distance, BAz etc.)
"""
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from obspy.imaging.beachball import beach
from obspy.geodetics import gps2dist_azimuth

sys.path.append('../modules')
from getdata import pathnames
from getdata import get_moment_tensor
from procmod import myround

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.linewidth'] = 2

# ================================= UTILITIES ================================= 
def build_colormap(array):
    """build a custom range colormap 
    """
    # build colormap for misfits
    vmax = myround(np.nanmax(array),base=1,choice='up')
    vmin = myround(np.nanmin(array),base=1,choice='down')
    # vmax = np.nanmax(array)
    # vmin = np.nanmin(array)
    norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    cmap = cm.plasma
    colormap = cm.ScalarMappable(norm=norm,cmap=cmap)
    
    return colormap

# ============================ DRAWER FUNCTIONS ================================
def trace_trench(m):
    """trace the hikurangi trench on a basemap object 'm'
    """
    trenchcoordspath = pathnames()['data'] + \
                                        'FAULTS/hikurangiTrenchCoords.npz'
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

def onshore_offshore_faults(m):
    """plot onshore and offshore fault coordinate files
    """
    onshore_fault_path = pathnames()['data'] + "FAULTS/onshoreFaultCoords.npz"
    offshore_fault_path = pathnames()['data'] + "FAULTS/offshoreFaultCoords.npz"

    onshore = np.load(onshore_fault_path)
    offshore = np.load(offshore_fault_path)

    for shore in [onshore,offshore]:
        lats = shore['LAT']
        lons = shore['LON']
        faults = shore['FAULT']

        for i in range(faults.min(),faults.max()+1,1):
            indices = np.where(faults==i)
            x,y = m(lons[indices],lats[indices])
            m.plot(x,y,'--',linewidth=0.5,color='k',zorder=2,alpha=0.25)

def event_beachball(m,MT):
    """plot event beachball on basemap 'm' object for a given geonet event_id
    """
    try:
        eventx,eventy = m(MT['Longitude'],MT['Latitude'])
        FM = [MT['strike2'],MT['dip2'],MT['rake2']]

        b = beach(FM,xy=(eventx,eventy),
                     width=2.5E4,
                     linewidth=1,
                     facecolor='r')
        b.set_zorder(1000)
        ax = plt.gca()
        ax.add_collection(b)
        return True
    except Exception as e:
        print('No moment tensor information found, beachball not available')
        plt.plot()
        return False
    

# ============================== MAP INITIATION ================================
def initiate_basemap(map_corners=[-50,-32.5,165,180]):
    """set up local map of NZ to be filled
    default map corners give a rough box around new zealand
    etopo_else_flat_map: bool to control continent fill with either low res topo
    or solid line representation
    """
    # flat map configuration
    continent_color = 'w'
    lake_color = 'w'

    # initiate map
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
    m.drawparallels(np.arange(int(map_corners[0]),int(map_corners[1]),1),
                                    labels=[1,0,0,0], linewidth=0.0)
    m.drawmeridians(np.arange(int(map_corners[2]),int(map_corners[3])+1,1),
                                    labels=[0,0,0,1], linewidth=0.0)

    return m

# ============================== MAP POPULATION ================================
def populate_basemap(m,lats,lons,names=None):
    """fill map with latitude/longitude pairs, i.e. stations, events
    """
    if names is None: names = []

    X,Y = m(lons,lats)
    scatter = m.scatter(X,Y,
                        marker='v',
                        color='w',
                        edgecolor='k',
                        s=60,
                        zorder=5)

    if len(names) != 0:
        for n_,x_,y_ in zip(names,X,Y):
            plt.annotate(n_,xy=(x_,y_),xytext=(x_,y_),zorder=6,fontsize=8.5)
            
def event_info_anno(m,srcrcvdict):
    """annotate event information into hard coded map area
    """
    # change annotate template depending on what type of map is produced
    if srcrcvdict['backazimuth']:
        annotemplate = ("{id} / {sta}\n{date}\n{type}={mag:.2f}"
            "\nDepth(km)={depth:.2f}\nDist(km)={dist:.2f}\nBAz(deg)={baz:.2f}")
    else:    
        annotemplate = ("{id} / {sta}\n{date}\n"
                        "{type}={mag:.2f}\nDepth(km)={depth:.2f}")
    
    plt.annotate(s=annotemplate.format(id=srcrcvdict['event_id'],
                                       sta=srcrcvdict['station'],
                                       date=srcrcvdict['date'],
                                       stalat=srcrcvdict['sta_lat'],
                                       stalon=srcrcvdict['sta_lon'],
                                       evlat=srcrcvdict['ev_lat'],
                                       evlon=srcrcvdict['ev_lon'],
                                       type=srcrcvdict['magnitude_type'],
                                       mag=srcrcvdict['magnitude'],
                                       depth=srcrcvdict['depth'],
                                       dist=srcrcvdict['distance'],
                                       baz=srcrcvdict['backazimuth']),
                 xy=((m.xmax-m.xmin)*0.725,(m.ymax-m.ymin)*0.0375),
                 # weight='bold',
                 multialignment='right',
                 fontsize=10)

def source_receiver(m,event,inv=None):
    """determine source receiver parameters such as great circle distance,
    backazimuth etc., plot the source and receiver with a line. have the ability
    to not include an inventory, which will just plot the event then
    source = A, receiver = B
    """
    # get event information
    event_id = event.resource_id.id.split('/')[1]
    MT = get_moment_tensor(event_id=event_id)
    event_lat,event_lon = (event.origins[0].latitude,
                           event.origins[0].longitude)
    # modify entries to dictionary
    origintime = event.origins[0].time
    origintime.precision = 0
    depth = event.origins[0].depth*1E-3
    # get magnitude M
    for magni in event.magnitudes:
        if magni.magnitude_type == "M":
            magnitude = magni.mag
            magnitude_type = magni.magnitude_type
        
    # if you want to plot a station as well
    if inv:
        station = inv[0][0].code
        sta_lat,sta_lon = inv[0][0].latitude,inv[0][0].longitude
        GCDist,Az,BAz = gps2dist_azimuth(event_lat,event_lon,sta_lat,sta_lon)
        GCDist *= 1E-3
    else:
        station,sta_lat,sta_lon = "Misfits",None,None
        GCDist,Az,BAz = None,None,None

    # dictionary output for use in annotations
    srcrcvdict = {"station":station,
                  "sta_lat":sta_lat,
                  "sta_lon":sta_lon,
                  "ev_lat":event_lat,
                  "ev_lon":event_lon,
                  "event_id":event_id,
                  "distance":GCDist,
                  "backazimuth":BAz,
                  "date":origintime,
                  "depth":depth,
                  "magnitude":magnitude,
                  "magnitude_type":magnitude_type
                  }

    # connect source receiever with line and color receiver, plot event
    event_x,event_y = m(event_lon,event_lat)
    if inv:
        sta_x,sta_y = m(sta_lon,sta_lat)
        X,Y = [event_x,sta_x],[event_y,sta_y]
        m.plot(X,Y,'--',linewidth=1.1,c='k',zorder=2)
    else:
        X,Y = event_x,event_y
    m.scatter(X,Y,marker='v',color='r',edgecolor='k',s=80,zorder=6)

    beachballcheck = event_beachball(m,MT)
    if not beachballcheck:
        m.scatter(event_x,event_y,marker='o',color='r',edgecolor='k',
                                                        s=200,zorder=6)

    return srcrcvdict

def plot_misfits(f,m,cork,comp="Z",add_colorbar=True):
    """take a misfit dictionary from a Cork object and plot onto the basemap
    """

    stations,latitudes,longitudes = [],[],[]
    for sta in cork.stations:
        coordict = cork.ds.waveforms[sta].coordinates
        stations.append(sta)
        latitudes.append(coordict["latitude"])
        longitudes.append(coordict["longitude"])
    
    # break out misfit measures from cork object, build a colormap
    misfits = np.fromiter(cork.misfit_values.values(),dtype="float")
    colormap = build_colormap(misfits)
    
    # plot onto basemap object
    xs,ys = m(longitudes,latitudes)
    colors = []
    for sta,X,Y in zip(stations,xs,ys):
        stationcode = "{s}.HH{c}".format(s=sta,c=comp)
        try:
            misfit = cork.misfit_values[stationcode]
        except KeyError:
            continue
        colors.append(colormap.to_rgba(misfit))
    
    m.scatter(xs,ys,marker='v',color=colors,edgecolor='k',s=150,zorder=6)
    
    if add_colorbar:
        comp_dict = {"Z":"vertical","N":"north","E":"east",
                     "T":"transvserse","R":"radial"}
        colormap.set_array(misfits)
        cbar = f.colorbar(colormap,fraction=0.046,pad=0.04)
        cbar.set_label("{} misfit".format(comp_dict[comp]),
                                    rotation=270,labelpad=17,fontsize=14.5)
        

# ============================== MAP GENERATION ================================
def generate_map(event,inv,corners=[-42.5007,-36.9488,172.9998,179.5077],
                                                        faults=False,**kwargs):
    """initiate and populate a basemap object for New Zealands north island.
    Functionality to manually ignore stations based on user quality control
    Takes station coordinates and coloring from npz files
    Choice to annotate two stations which correspond to waveforms
    Calls beachball and trench tracer to populate basemap
    :type corners: list of floats
    :param corners: values for map corners to set bounds
     e.g. [lat_bot,lat_top,lon_left,lon_right]
    """
    m = initiate_basemap(map_corners=corners)
    srcrcvdict = source_receiver(m,event,inv)
    event_info_anno(m,srcrcvdict)
    if faults:
        trace_trench(m)
        onshore_offshore_faults(m)

    stationfile = pathnames()['data'] + 'STATIONXML/GEONET_AND_FATHOM.npz'
    stationlist = np.load(stationfile)

    populate_basemap(m,stationlist['LAT'],stationlist['LON'])
    scalelon,scalelat = 178.75,-37.2
    m.drawmapscale(scalelon,scalelat,scalelon,scalelat,100,
                                                yoffset=0.01*(m.ymax-m.ymin))

    return m
    
def generate_misfit_map(event_id,corners=[-42.5007,-36.9488,172.9998,179.5077],
                                                        faults=False,**kwargs):
    """take an event id associated with a PyASDF dataset that has been filled 
    out by the adjointBuilder and plot the misfits recorded at each station
    """
    import pyasdf
    from corkBoard import Cork
    
    f = plt.figure(figsize=(10,9.4),dpi=100)

    # use corkBoard to interact with pyASDF dataformat
    mycork = Cork(event_id)
    mycork.populate()
    mycork.get_srcrcv_information()
    mycork.collect_misfits()
        
    m = initiate_basemap(map_corners=corners)
    srcrcvdict = source_receiver(m,mycork.ds.events[0],inv=None)
    event_info_anno(m,srcrcvdict)
    if faults:
        trace_trench(m)
        onshore_offshore_faults(m)

    stationfile = pathnames()['data'] + 'STATIONXML/GEONET_AND_FATHOM.npz'
    stationlist = np.load(stationfile)

    populate_basemap(m,stationlist['LAT'],stationlist['LON'])
    plot_misfits(f,m,mycork)
    
    scalelon,scalelat = 178.75,-37.2
    m.drawmapscale(scalelon,scalelat,scalelon,scalelat,100,
                                                yoffset=0.01*(m.ymax-m.ymin))
                                                
    plt.show()

    return m

def _test_generate_map():
    """test map making functions with example data
    """
    from obspy import read_events, read_inventory

    # eventpath = pathnames()['data'] + "WINDOWTESTING/testevent.xml"
    invpath = pathnames()['data'] + "WINDOWTESTING/testinv.xml"

    # cat = read_events(eventpath)
    # event = cat[0]
    inv = read_inventory(invpath)
    
    m = generate_map_NOEVENT(inv,show=True,faults=True)
    # f,m = generate_map(event,inv,show=True)


if __name__ == "__main__":
    # _test_generate_map()
    generate_misfit_map("2014p240655")
