"""
Artificially scale the components of a moment tensor so that they sum up to
a desired Mw. Used for scaling tectonic events to the size of a declared test
so that it can be used in a simulation run
"""
from math import sqrt, log10
from obspy.imaging.beachball import beachball

PLOT = False
# Desired Mw
mw_desired = 4.0

# Original moment tensor values
mt_dict = dict(m_rr = 5.395e+23,
               m_tt = 5.395e+23,
               m_pp = 5.395e+23,
               m_rt = 0,
               m_rp = 0,
               m_tp = 0 
               )

# Check what the current Mw is
mt = list(mt_dict.values())
sum_mt_squared = sum([_ ** 2 for _ in mt])
m0 = (1 / sqrt(2)) * (sqrt(sum_mt_squared))
mw = (2/3) * log10(m0) - 10.7
print(f"current mw = {mw:.2f}")

# Calculate the scale factor to achieve the desired Mw
scale_factor = 10 ** ((3/2) * (mw_desired - mw))
print(f"calculated scale_factor = {scale_factor:.2f}")

# Apply the scale factor
mt_new = [_ * scale_factor for _ in mt]
sum_mt_squared_new = sum([_ ** 2 for _ in mt_new])
m0 = (1 / sqrt(2)) * (sqrt(sum_mt_squared_new))
mw_scaled = (2/3) * log10(m0) - 10.7
print(f"scaled mw = {mw_scaled:.2f}")

for i, (key, val) in enumerate(mt_dict.items()):
    suffix = key.split("_")[1]
    print(f"M{suffix}:      {val * scale_factor:13.6e}")

if PLOT:
    beachball(mt, size=100, linewidth=2, facecolor="b")
    beachball(mt_new, size=100, linewidth=2, facecolor="r")

