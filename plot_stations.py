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

vic_or_gns = 'GNS'

# choices: s (accelerometers) / z (seismometer) / d (durations)
choice = sys.argv[1].upper()
if choice not in ['S','Z','D']:
    sys.exit("\tChoice should be s(acc), z(seis) or d(dur)")
tmin = 1
tmax = 50
perc = 90
sta_dict = {'Z':'*Z','S':'*S','D':'*S'}
chan_dict = {'Z':'HH*','S':'BN*','D':'BN*'}
color_dict = {'Z':'b','S':'r','D':'r'}
text_dict = {'Z':'Seismometer','S':'Accelerometer','D':'Accelerometer'}
label_dict = {'Z':True,'S':True,'D':False}


c = Client('GEONET')
# 2 degree radius around Te Urewera
north_island = [-45,-35,175,180]
new_zealand = [-50,-35,165,180]
inv = c.get_stations(network='NZ',
                    station=sta_dict[choice],
                    minlatitude=-41,
                    maxlatitude=-39,
                    minlongitude=175,
                    maxlongitude=180,
                    level="station")

fig = inv.plot(label=label_dict[choice],
        marker='.',
        projection="local",
        resolution = "i",
        color=color_dict[choice],
        size=6,
        show=False)

# ========= plot RDF station array =========
stations,lats,lons = [],[],[]
with open(pathnames(vic_or_gns)['rdf'] + 'rdf_locations.txt','r') as f:
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
text_file = (pathnames(vic_or_gns)['spectral'] +
                                'duration/{t0}-{t1}s_amp.txt'.format(
                                                                t0=tmin,
                                                                t1=tmax))
stations,durations,lats,lons = [],[],[],[]
with open(text_file,'r') as f:
    for line in f:
        if "ID" in line:
            station_temp = line.split('.')[1]
            stations.append(station_temp)
        elif "Duration" in line:
            line_temp = line.strip()
            duration_temp = float(line.split(' ')[1])
            durations.append(duration_temp)

# grab latitude and longitude information
for sta in stations:
    sta_info = inv.select(station=sta)[0][0]
    lats.append(sta_info.latitude)
    lons.append(sta_info.longitude)

# set colors
D = [(_/max(durations))**2 for _ in durations]
cmap = plt.get_cmap('YlGn')
colors = cmap(D)

# plot onto basemap
x,y = fig.bmap(lons,lats)
scatter = fig.bmap.scatter(x,y,marker='o',s=100,zorder=1000,c=colors)
for i,(sta,dur) in enumerate(zip(stations,durations)):
    plt.annotate("{}/{}s".format(sta,dur),
                    xy=(x[i],y[i]),
                    xytext=(x[i]+5000,y[i]+5000),
                    fontsize=8)

plt.title('Geonet {name} w/ shaking durations | Filt: {Fl}-{Fh} s'.format(
            name=text_dict[choice],
            Fl=tmin,
            Fh=tmax))

plt.show()
