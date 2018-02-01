"""use obspy to create station maps
"""
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from obspy.clients.fdsn import Client

# ignore warnings
import warnings
warnings.filterwarnings("ignore",category=mpl.cbook.mplDeprecation)

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
lat = -38.468
lon = 177.217849

c = Client('GEONET')
# 2 degree radius around Te Urewera
inv = c.get_stations(network='NZ',
                    station=sta_dict[choice],
                    channel=chan_dict[choice],
                    starttime='2015-01-01',
                    endtime='2016-01-01',
                    latitude=lat,
                    longitude=lon,
                    maxradius=5)

fig = inv.plot(label=label_dict[choice],
        marker='.',
        projection="local",
        resolution = "i",
        color=color_dict[choice],
        size=6,
        show=False)

if choice != 'D':
    plt.show()
    sys.exit()

# manual plot scatterpoints - really hacked together
text_file = '/seis/prj/fwi/bchow/spectral/duration/{t0}-{t1}s_amp.txt'.format(
                                                                t0=tmin,
                                                                t1=tmax)
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
