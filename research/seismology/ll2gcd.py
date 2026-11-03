"""
Lat Lon to Great Circle Distance
"""
import sys
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import kilometer2degrees

try:
    lat1, lon1, lat2, lon2 = sys.argv[1:]
except Exception as e:
    print("ll2gcd <lat1> <lon1> <lat2> <lon2>")

gcd, az, baz = gps2dist_azimuth(float(lat1), float(lon1),
                                float(lat2), float(lon2))
gcd_deg = kilometer2degrees(gcd)

print(f"GCD     = {gcd*1E-3:.2f} km")
print(f"GCD     = {gcd_deg*1E-3:.2f} deg")
print(f"AZ A->B = {az:.2f} deg")
print(f"AZ B->A = {baz:.2f} deg")
                                
