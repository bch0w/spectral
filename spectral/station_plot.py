"""use obspy to create station maps
"""
import sys
sys.path.append('../modules/')
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from getdata import pathnames, get_moment_tensor
from obspy import read_inventory, read_events
from matplotlib.patches import Polygon
from obspy.geodetics import locations2degrees
from obspy.imaging.beachball import beach

# ignore warnings
import warnings
warnings.filterwarnings("ignore",category=mpl.cbook.mplDeprecation)

def event_beachball(eventid,fig):
    """plt event beachball on figure object
    """
    MT = get_moment_tensor(eventid)
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

def rdf_plot(fig):
    """plot rdf on map - not tested
    """
    stations,lats,lons = [],[],[]
    with open(pathnames()['rdf'] + 'rdf_locations.txt','r') as f:
        station_lines = f.readlines()

    for lines in station_lines[1:]:
        splitlines = lines.split(',')
        # stations.append("{n}/{s}".format(n=splitlines[0],s=splitlines[1]))
        stations.append(splitlines[0])
        lats.append(float(splitlines[2]))
        lons.append(float(splitlines[3]))

    # plot onto basemap
    x,y = fig.bmap(lons,lats)
    scatter = fig.bmap.scatter(x,y,marker='o',s=40,zorder=1000,c='g')
    for i,sta in enumerate(stations):
        plt.annotate(sta,xy=(x[i],y[i]),
                        xytext=(x[i]+2500,y[i]+2500),
                        fontsize=9,
                        fontweight='bold')
    import ipdb;ipdb.set_trace()

    map_limx,map_limy = fig.bmap([175.75,178],[-39,-40.75])
    plt.xlim(map_limx)
    plt.ylim(map_limy)

def low_velocity_wedge(fig):
    """plot low velocity wedge on map
    """
    nw_lat,nw_lon = -39.8764, 178.5385
    ne_lat,ne_lon = -39.0921, 177.1270
    se_lat,se_lon = -37.5630, 178.5248
    sw_lat,sw_lon = -38.3320, 179.9157
    x1,y1 = fig.bmap(nw_lon,nw_lat)
    x2,y2 = fig.bmap(ne_lon,ne_lat)
    x3,y3 = fig.bmap(se_lon,se_lat)
    x4,y4 = fig.bmap(sw_lon,sw_lat)
    poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],
                    facecolor='red',edgecolor='k',
                    linewidth=1,alpha=0.2)
    ax = plt.gca()
    ax.add_patch(poly)

# ================================== MAIN ====================================
# command line choices
choice = sys.argv[1].upper() #s,z,d
if choice == 'D':
    event_id = sys.argv[2] #2014p240655,2015p822263,2016p842451,2922302
    component = sys.argv[3].lower() #vertical,horizontal,groundmotion
if choice not in ['S','Z','D']:
    sys.exit("\tChoice should be s(acc), z(seis) or d(dur)")

# dictionaries for dynamic choices
comp_dic = {"vertical":1,
            "horizontal":2,
            "groundmotion":3}
sta_dict = {'Z':'*Z','S':'*S','D':'*Z'}
chan_dict = {'Z':'HH*','S':'*N*','D':'HH*'}
color_dict = {'Z':'b','S':'r','D':'r'}
text_dict = {'Z':'Seismometer','S':'Accelerometer','D':'Seismometer'}
label_dict = {'Z':True,'S':True,'D':False}

# grab station information and plot
inv = read_inventory(pathnames()["data"] +
                            "STATIONXML/new_zealand_geonet_stations.xml")
inv += read_inventory(pathnames()["data"] +
                            "STATIONXML/hobitss_stations.xml")

fig = inv.plot(marker='.',
        projection="local",
        resolution = "i",
        color=color_dict[choice],
        size=25,
        show=False,
        continent_fill_color="white",
        water_fill_color="white",
        color_per_network=True,
        label=False)
        # label=label_dict[choice])

low_velocity_wedge(fig)
for EID in ["2014p864702","2014p715167","2014p240655",
                                        "2015p822263","2014p051675"]:
    event_beachball(EID,fig)

if choice != 'D':
    plt.show()
    sys.exit()

# ==================== PROCESSING FOR MAP =====================================
event_pathname = pathnames()["data"] + event_id
event = read_events(event_pathname,format="QUAKEML")[0]
event_lat,event_lon = event.origins[0].latitude,event.origins[0].longitude
event_mag = round(event.magnitudes[0].mag,2)

# retrieve values from text file
text_file = pathnames()['spectral']+'duration/{}_5-30_amplitudes.txt'.format(
                                                                    event_id)
durations,durstations,lats,lons = [],[],[],[]
with open(text_file,'r') as f:
    for line in f:
        eachline = line.split(' ')
        if "#" in eachline[0]:
            continue
        durstations.append(eachline[0])
        durations.append(float(eachline[comp_dic[component]]))

# grab latitude and longitude information
for sta in durstations:
    sta_info = inv.select(station=sta)[0][0]
    lats.append(sta_info.latitude)
    lons.append(sta_info.longitude)

# set colors
D = [(_/max(durations)) for _ in durations]
cmap = plt.get_cmap('plasma')
colors = cmap(D)

# determine station distances as colors
distances = []
for La,Lo in zip(lats,lons):
    dist = locations2degrees(La,Lo,event_lat,event_lon)
    dist *= (2.0 * 6371 * np.pi / 360.0)
    distances.append(dist)
dist_plot = distances
distances = [_/max(distances) for _ in distances]
cmap2 = plt.get_cmap('plasma_r')
distance_colors = cmap2(distances)

# ================================= PLOT ON BASEMAP ===========================
# plot stations
map_x,map_y = fig.bmap(lons,lats)
scatter = fig.bmap.scatter(map_x,map_y,
                            marker='v',
                            s=150,
                            zorder=1000,
                            c=colors,
                            edgecolor='k')

# annotate stations
for i,(sta,dur,ds) in enumerate(zip(durstations,durations,dist_plot)):
    plt.annotate("{0}({1}mm/s)[{2}km]".format(sta,round(dur,2),int(ds)),
                    xy=(map_x[i],map_y[i]),
                    xytext=(map_x[i]+5000,map_y[i]+5000),
                    fontsize=10,
                    weight='bold',
                    zorder=500)



# connect event and stations with colored lines
for X,Y,DC in zip(map_x,map_y,distance_colors):
    dist = locations2degrees(La,Lo,event_lat,event_lon)
    dist *= (2.0 * 6371 * np.pi / 360.0)
    plt.annotate('', xy=(X,Y),
                    xycoords='data',
                    xytext=(eventx,eventy),
                    textcoords='data',
                    arrowprops=dict(arrowstyle="-",color=DC),
                    alpha=0.25,
                    zorder=1)

# set limits via lat lon coordinates
LLC = fig.bmap(173,-41)
URC = fig.bmap(180,-37)
plt.xlim([LLC[0],URC[0]])
plt.ylim([LLC[1],URC[1]])

# final figure adjustments
plt.title('{name} {comp} peak amplitudes | {Eid}'.format(
            comp=component,
            name=text_dict[choice],
            Eid=event_id))
fname = "{comp}_{Eid}.png".format(comp=component,Eid=event_id)
fullname = pathnames()['plots'] + "durations/{}".format(fname)

plt.savefig(fullname,dpi=200)

plt.show()
