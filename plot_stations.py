"""use obspy to create station maps
"""
import sys
import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client

choice = sys.argv[1].upper()
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

fig = inv.plot(label=True,
        marker='^',
        projection="local",
        resolution = "i",
        color=color_dict[choice],
        size=6,
        show=False)

plt.title('Geonet {name} - 2 deg radius around ({lat},{lon})'.format(
            name=text_dict[choice],
            lat=lat,
            lon=lon))
plt.show()
