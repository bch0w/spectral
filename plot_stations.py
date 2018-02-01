"""use obspy to create station maps
"""
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from obspy.clients.fdsn import Client

# ignore warnings
import warnings
warnings.filterwarnings("ignore",category=mpl.cbook.mplDeprecation)

choice = sys.argv[1].upper()
tmin = 1
tmax = 50
perc = 90
chan_dict = {'Z':'HH*','S':'BN*'}
color_dict = {'Z':'b','S':'r'}
text_dict = {'Z':'Seismometer','S':'Accelerometer'}
lat = -38.468
lon = 177.217849

c = Client('GEONET')
# 2 degree radius around Te Urewera
inv = c.get_stations(network='NZ',
                    station='*{}'.format(choice),
                    channel=chan_dict[choice],
                    starttime='2015-01-01',
                    endtime='2016-01-01',
                    latitude = lat,
                    longitude = lon,
                    maxradius=2)

fig = inv.plot(label=False,
        marker='.',
        projection="local",
        resolution = "i",
        color=color_dict[choice],
        size=6,
        show=False)

# manual plot scatterpoints - really hacked together
text_file = '/seis/prj/fwi/bchow/spectral/durations/{t0}-{t1}s_{p}p.txt'.format(
                                                                t0=tmin,
                                                                t1=tmax,
                                                                p=perc)
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
