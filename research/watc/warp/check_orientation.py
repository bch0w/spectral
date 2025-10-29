"""
Look at nearby earthquakes to get orientation, very brute force
"""
import sys
import matplotlib.pyplot as plt
import numpy as np
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth

choice = 12013854

if choice == 12013854:
    # https://ds.iris.edu/wilber3/find_stations/12013854
    start = UTCDateTime("2025-08-31T16:05:26")
    end = UTCDateTime("2025-08-31T16:05:30")
    ev_lat = 57.459
    ev_lon = -153.3107
elif choice == 12014072:
    # https://ds.iris.edu/wilber3/find_stations/12014072
    start = UTCDateTime("2025-09-01T11:02:30Z")
    end = UTCDateTime("2025-09-01T11:02:35Z")
    ev_lat = 56.2794
    ev_lon = -151.9837

c = Client("IRIS")
code = "II.KDAK.10.BH?"
net, sta, loc, cha = code.split(".")  # 00=Borehole; 10=surface
st = c.get_waveforms(network=net, station=sta, location=loc, channel=cha,
                     starttime=start, endtime=end)
inv = c.get_stations(network=net, station=sta, location=loc, channel=cha,
                     starttime=start, endtime=end, level="response")

sta_lat = inv[0][0].latitude
sta_lon = inv[0][0].longitude

dist, az_exp, baz_exp = gps2dist_azimuth(ev_lat, ev_lon, sta_lat, sta_lon)

# BH1 azimuth = 155.20* clockwise from North
# BH2 azimuth = 245.30* clockwise from North
st.rotate(method="->ZNE", inventory=inv)

z, r, t = [], [], []
backazimuths = np.arange(0, 360, 1)
for baz in backazimuths:
    st_r = st.copy()
    st_r.rotate(method="NE->RT", back_azimuth=baz)
    r.append(np.amax(st_r.select(component="R")[0].data))
    t.append(np.amax(st_r.select(component="T")[0].data))
    z.append(np.amax(st_r.select(component="Z")[0].data))

r_max = backazimuths[np.where(r==max(r))][0]
t_max = backazimuths[np.where(t==max(t))][0]
z_max = backazimuths[np.where(z==max(z))][0]

plt.plot(backazimuths, r, "ro-", label=f"radial ({r_max:.2f})")
plt.plot(backazimuths, t, "bo-", label=f"transverse ({t_max:.2f})")
# plt.plot(backazimuths, z, "ko-", label=f"vertical ({z_max:.2f})")
plt.title(f"{code} {start}\n"
          f"Expected BAz = {baz_exp:.2f}; "
          f"delta_BAz = {r_max - baz_exp:.2f}")
plt.grid()
plt.legend()
plt.xlabel("BAz")
plt.ylabel("Max P-wave Amp. [counts]")
plt.savefig(f"{code}_{choice}_orientation_check.png")

plt.show()


