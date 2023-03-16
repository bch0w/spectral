"""
Rotate T component adjoint sources, which would have been created from the 
comparison of TT SGF and EGF, back to E and N components for use in adjoint
simulations.

Note that this will generate FOUR sets of adjoint sources, E_E, N_E, E_N, N_N,
where the subscript corresponds to the forward simulation force component. 
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
for net in inv:
    for sta in net:
        n = net.code
        s = sta.code

        # Time and amplitude for adjoint source (u_tt)
        try:
            t_arr, u_tt = np.loadtxt(f"T/{n}.{s}.BXT.adj", usecols=[0, 1]).T
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

        # Rotate T into various E and N combinations
        u_ee = +1 * np.cos(theta) * np.cos(theta_p) * u_tt
        u_en = -1 * np.cos(theta) * np.sin(theta_p) * u_tt 
        u_ne = -1 * np.sin(theta) * np.cos(theta_p) * u_tt 
        u_nn = +1 * np.sin(theta) * np.sin(theta_p) * u_tt 


        np.savetxt(f"EE/{n}.{s}.BXE.adj", np.vstack((t_arr, u_ee)).T, 
                   fmt="%11.6f%21.7E")
        np.savetxt(f"EN/{n}.{s}.BXE.adj", np.vstack((t_arr, u_en)).T, 
                   fmt="%11.6f%21.7E")
        np.savetxt(f"NE/{n}.{s}.BXN.adj", np.vstack((t_arr, u_ne)).T, 
                   fmt="%11.6f%21.7E")
        np.savetxt(f"NN/{n}.{s}.BXN.adj", np.vstack((t_arr, u_nn)).T, 
                   fmt="%11.6f%21.7E")

