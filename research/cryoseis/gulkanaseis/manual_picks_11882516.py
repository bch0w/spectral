import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime as UTC
from obspy.geodetics import gps2dist_azimuth

# Manually picked P-wave arrival times
p_picks = {
        "100": UTC("2024-09-09T01:40:52.37586"),
        "101": UTC("2024-09-09T01:40:52.37598"),
        "102": UTC("2024-09-09T01:40:52.37599"),
        "103": UTC("2024-09-09T01:40:52.37599"),
        "104": UTC("2024-09-09T01:40:52.37600"),
        "105": UTC("2024-09-09T01:40:52.37999"),
        "106": UTC("2024-09-09T01:40:52.38402"),
        "107": UTC("2024-09-09T01:40:52.38042"),
        "108": UTC("2024-09-09T01:40:52.37994"),
        "109": UTC("2024-09-09T01:40:52.39598"),
        "110": UTC("2024-09-09T01:40:52.38794"),
        "111": UTC("2024-09-09T01:40:52.39200"),
        }

# Manually picked S-wave arrival times
s_picks = {
        "100": UTC("2024-09-09T01:41:11.89588"),
        "101": UTC("2024-09-09T01:41:11.89603"),
        "102": UTC("2024-09-09T01:41:11.89587"),
        "103": UTC("2024-09-09T01:41:11.89993"),
        "104": UTC("2024-09-09T01:41:11.90002"),
        "105": UTC("2024-09-09T01:41:11.89596"),
        "106": UTC("2024-09-09T01:41:11.90796"),
        "107": UTC("2024-09-09T01:41:11.90780"),
        "108": UTC("2024-09-09T01:41:11.89977"),
        "109": UTC("2024-09-09T01:41:11.91195"),
        "110": UTC("2024-09-09T01:41:11.91193"),
        "111": UTC("2024-09-09T01:41:11.90791"),
        }

# GPS locations from Garmin
locations = {
    "100": (63.265612, -145.415025),
    "101": (63.265362, -145.41407),
    "102": (63.265448, -145.41381),
    "103": (63.265203, -145.415355),
    "104": (63.265062, -145.414323),
    "105": (63.26511, -145.413887),
    "106": (63.26495, -145.415472),
    "107": (63.264683, -145.414657),
    "108": (63.264878, -145.414107),
    "109": (63.26392, -145.41601),
    "110": (63.26466, -145.41557),
    "111": (63.264278, -145.415858),
}

stations = p_picks.keys()

# Set reference, we know that 100 gets the signal first
ref_station = "100"
ref_p_pick = p_picks[ref_station]
ref_s_pick = s_picks[ref_station]
ref_lat, ref_lon = locations[ref_station]

# Loop through stations and calculate
dists_m, diffs_p, diffs_s, p_minus_s = [], [], [], []
apparent_p, apparent_s = [], []

for sta in stations:
    # Calculate station distances
    lat, lon = locations[sta]
    dist_m, *_ = gps2dist_azimuth(ref_lat, ref_lon, lat, lon)
    dists_m.append(dist_m)

    # Get arrival time differences across array
    diff_p = p_picks[sta] - ref_p_pick
    diffs_p.append(diff_p)

    diff_s = s_picks[sta] - ref_s_pick
    diffs_s.append(diff_s)

    # Calculate velocity based on traveltime across array
    if sta != ref_station:
        p_vel = np.abs(dist_m / diff_p) * 1E-3  # km/s
        s_vel = np.abs(dist_m / diff_s) * 1E-3  # km/s
    else:
        p_vel, s_vel = 0, 0
    apparent_p.append(p_vel)
    apparent_s.append(s_vel)

    # Get P-S arrival time difference
    p_minus_s.append(s_picks[sta] - p_picks[sta])

breakpoint()



