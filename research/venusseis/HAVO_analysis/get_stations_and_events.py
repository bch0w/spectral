import matplotlib.pyplot as plt
from obspy.clients.fdsn import Client
from obspy import UTCDateTime


# Hawaiin Archipelago
lon_min = -162.8
lon_max = -150.0
lat_min = 15.9
lat_max = 24.4

mag_min = 2.0
mag_max = 6.0

# Last decades of seismicity
starttime = UTCDateTime("2014-01-01T00:00:00")
endtime = UTCDateTime("2024-10-27T00:00:00")

# Get stations and data
c = Client("IRIS")
inv = c.get_stations(starttime=starttime, endtime=endtime, network="*", 
                     channel="HH?")
cat = c.get_events(starttime=starttime, endtime=endtime, 
                   minlongitude=lon_min, maxlongitude=lon_max,
                   minlatitude=lon_min, maxlatitude=lat_max, 
                   minmagnitude=mag_min, maxmagnitude=mag_max,
                   )

print(f"{len(inv)} stations")
print(f"{len(cat)} events")


f, ax = plt.subplot()
cat.plot(projection="local", resolution="c", show=False, fig=f)
inv.plot(fig=f)
plt.show()
