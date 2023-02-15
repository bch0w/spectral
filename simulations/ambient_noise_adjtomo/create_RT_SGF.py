"""
From SPECFEM outputs, calculate TT SGF from NN + EE SGF following equations 
from Wang et al. (2019)
"""
import numpy as np
from glob import glob
from obspy.geodetics import gps2dist_azimuth
from pysep.utils.io import read_stations


# Set parameters here 
stations_file = "../DATA/STATIONS"

# Used to grab station coordinates
inv = read_stations(stations_file)

# Find information for master station
# Hardcode location of: IU.COLA
lat_a = 64.8736
lon_a = -147.8616

if not lat_a:
    station_a = inv.select(network="IU", station="COLA")
    lat_a = station_a[0][0].latitude
    lon_a = station_a[0][0].longitude

# Dummy station to grab time axis for. Assumed all will have same t-axis
t_arr = np.loadtxt(glob(f"N/*.ascii")[0], usecols=0)
for net in inv:
    for sta in net:
        n = net.code
        s = sta.code

        # Some combinations don't exist. Ignore these
        # 1st letter = point source direction; 2nd letter = recorded component
        try:
            u_nn = np.loadtxt(f"N/{n}.{s}.BXN.sem.ascii", usecols=1)
            u_ne = np.loadtxt(f"N/{n}.{s}.BXE.sem.ascii", usecols=1)
            u_ee = np.loadtxt(f"E/{n}.{s}.BXE.sem.ascii", usecols=1)
            u_en = np.loadtxt(f"E/{n}.{s}.BXN.sem.ascii", usecols=1)
        except FileNotFoundError:
            continue

        print(f"{n}.{s}")
        lat_b = sta.latitude
        lon_b = sta.longitude

        # See Fig 1 from Wang et al. (2019) for theta and theta' def.
        _, az, baz = gps2dist_azimuth(lat1=lat_a, lon1=lon_a, 
                                      lat2=lat_b, lon2=lon_b) 

        # theta != theta' for spherical Earth but will be very close
        theta = np.deg2rad(az)
        theta_p = np.deg2rad((baz - 180) % 360)  # theta' (p for prime)


        # Rotate N and E SGF to TT SGF
        u_tt = (+1 * np.cos(theta) * np.cos(theta_p) * u_ee 
                -1 * np.cos(theta) * np.sin(theta_p) * u_ne
                -1 * np.sin(theta) * np.cos(theta_p) * u_en
                +1 * np.sin(theta) * np.sin(theta_p) * u_nn)

        # Rotate N and E SGF to RR SGF (will not be used)
        u_rr = (+1 * np.sin(theta) * np.sin(theta_p) * u_ee 
                -1 * np.sin(theta) * np.cos(theta_p) * u_ne
                -1 * np.cos(theta) * np.sin(theta_p) * u_en
                +1 * np.cos(theta) * np.cos(theta_p) * u_nn)


        np.savetxt(f"T/{s}.{n}.BXT.sem.ascii", np.vstack((t_arr, u_tt)).T, 
                   fmt="%11.6f%21.7E")
        np.savetxt(f"R/{s}.{n}.BXR.sem.ascii", np.vstack((t_arr, u_rr)).T,
                   fmt="%11.6f%21.7E")
