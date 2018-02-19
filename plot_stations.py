"""use obspy to create station maps
"""
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from getdata import pathnames
from obspy.clients.fdsn import Client
from matplotlib.patches import Polygon

# ignore warnings
import warnings
warnings.filterwarnings("ignore",category=mpl.cbook.mplDeprecation)

# command line choices
choice = sys.argv[1].upper() #s,z,d
event_id = sys.argv[2] #2014p240655,2015p82263,2016p892721,2017p059122
component = sys.argv[3].lower() #vertical,horizontal,groundmotion
if choice not in ['S','Z','D']:
    sys.exit("\tChoice should be s(acc), z(seis) or d(dur)")

# dictionaries for dynamic choices
comp_dic = {"vertical":1,
            "horizontal":2,
            "groundmotion":3}
sta_dict = {'Z':'*Z','S':'*S','D':'*Z'}
chan_dict = {'Z':'*H*','S':'*N*','D':'HH*'}
color_dict = {'Z':'b','S':'r','D':'r'}
text_dict = {'Z':'Seismometer','S':'Accelerometer','D':'Seismometer'}
label_dict = {'Z':True,'S':True,'D':False}


# grab station information
c = Client('GEONET')
north_island = [-45,-35,175,180]
new_zealand = [-50,-35,165,180]
inv = c.get_stations(network='NZ',
                    station=sta_dict[choice],
                    minlatitude=-45,
                    maxlatitude=-35,
                    minlongitude=175,
                    maxlongitude=180,
                    level="station")

# plot stations
fig = inv.plot(label=label_dict[choice],
        marker='.',
        projection="local",
        resolution = "i",
        color=color_dict[choice],
        size=6,
        show=False)

# ========= plot RDF station array =========
    # stations,lats,lons = [],[],[]
    # with open(pathnames()['rdf'] + 'rdf_locations.txt','r') as f:
    #     station_lines = f.readlines()
    #
    # for lines in station_lines[1:]:
    #     splitlines = lines.split(',')
    #     # stations.append("{n}/{s}".format(n=splitlines[0],s=splitlines[1]))
    #     stations.append(splitlines[0])
    #     lats.append(float(splitlines[2]))
    #     lons.append(float(splitlines[3]))
    #
    # # plot onto basemap
    # x,y = fig.bmap(lons,lats)
    # scatter = fig.bmap.scatter(x,y,marker='o',s=40,zorder=1000,c='g')
    # for i,sta in enumerate(stations):
    #     plt.annotate(sta,xy=(x[i],y[i]),
    #                     xytext=(x[i]+2500,y[i]+2500),
    #                     fontsize=9,
    #                     fontweight='bold')

# manual set bounds
# map_limx,map_limy = fig.bmap([175.75,178],[-39,-40.75])
# plt.xlim(map_limx)
# plt.ylim(map_limy)

# ========= plot low velocity overlay with a polygon =========
nw_lat,nw_lon = -37.9281, 178.1972
ne_lat,ne_lon = -38.7006, 179.5930
se_lat,se_lon = -40.0229, 178.4043
sw_lat,sw_lon = -39.2371, 176.9909
x1,y1 = fig.bmap(nw_lon,nw_lat)
x2,y2 = fig.bmap(ne_lon,ne_lat)
x3,y3 = fig.bmap(se_lon,se_lat)
x4,y4 = fig.bmap(sw_lon,sw_lat)
poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],
                facecolor='red',edgecolor='k',
                linewidth=1,alpha=0.2)
plt.gca().add_patch(poly)

if choice != 'D':
    plt.show()
    sys.exit()

# ========= manual plot scatterpoints - really hacked together =========
text_file = pathnames()['spectral']+'duration/{}.txt'.format(event_id)
durations,durstations,lats,lons = [],[],[],[]
with open(text_file,'r') as f:
    for line in f:
        eachline = line.split(' ')
        if "time" in eachline[0]:
            station = eachline[0].split('_')[0]
            durstations.append(station)
            durations.append(float(eachline[comp_dic[component]]))
        else:
            continue

# grab latitude and longitude information
for sta in durstations:
    sta_info = inv.select(station=sta)[0][0]
    lats.append(sta_info.latitude)
    lons.append(sta_info.longitude)

# set colors
D = [(_/max(durations))**4 for _ in durations]
cmap = plt.get_cmap('YlGn')
colors = cmap(D)

# plot onto basemap
x,y = fig.bmap(lons,lats)
scatter = fig.bmap.scatter(x,y,marker='o',s=100,zorder=1000,c=colors)
for i,(sta,dur) in enumerate(zip(durstations,durations)):
    plt.annotate("{}/{}".format(sta,dur),
                    xy=(x[i],y[i]),
                    xytext=(x[i]+5000,y[i]+5000),
                    fontsize=8)
# set limits based off map - only works for predefined map, for more dynamic
# call use fig.bmap(lon,lat)
plt.xlim([370000,756892])
plt.ylim([248581,563109])
plt.title('{name} {comp} component durations | {Eid}'.format(
            comp=component,
            name=text_dict[choice],
            Eid=event_id))
fname = "{comp}_{Eid}.png".format(comp=component,Eid=event_id)
fullname = pathnames()['plots'] + "durations/{}".format(fname)
plt.savefig(fullname,dpi=200)

plt.show()
