"""
Testing out calculationg of TT SGF from NN + EE SGF following equations from 
Wang et al. (2019)
"""
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from pysep.utils.io import read_sem


# IU.COLA
lat_a = 64.8736
lon_a = -147.862

# TA.TOLK
lat_b = 68.6408
lon_b = -149.5724

# To see how theta and theta' are defined, see Fig 1 from Wang et al. (2019)
_, az, baz = gps2dist_azimuth(lat1=lat_a, lon1=lon_a, lat2=lat_b, lon2=lon_b) 

# theta != theta' for spherical Earth but will be very close
theta = np.deg2rad(az)
theta_p = np.deg2rad((baz - 180) % 360)  # theta' (p for prime)
print(f"theta={theta:.2f}; theta'={theta_p:.2f}")

# 1st letter = point source direction; 2nd letter = recorded component
# so e.g.: ne is North pointing point source recorded on the East component
time = np.loadtxt("NN_IU_COLA_TA_TOLK.sem.ascii", usecols=0)

u_nn = np.loadtxt("NN_IU_COLA_TA_TOLK.sem.ascii", usecols=1)
u_ne = np.loadtxt("NE_IU_COLA_TA_TOLK.sem.ascii", usecols=1)
u_ee = np.loadtxt("EE_IU_COLA_TA_TOLK.sem.ascii", usecols=1)
u_en = np.loadtxt("EN_IU_COLA_TA_TOLK.sem.ascii", usecols=1)

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


np.savetxt("TT_IU_COLA_TA_TOLK.sem.ascii", np.vstack((time, u_tt)).T)
