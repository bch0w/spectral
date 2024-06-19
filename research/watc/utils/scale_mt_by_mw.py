"""
Artificially scale the components of a moment tensor so that they sum up to
a desired Mw. Used for scaling tectonic events to the size of a declared test
so that it can be used in a simulation run
"""
from math import sqrt, log10
from obspy.imaging.beachball import beachball

# Desired Mw
mw_desired = 6.3

# Original moment tensor values
mt_dict = dict(m_rr =  2.280000e+23,
               m_tt =  1.550000e+24,
               m_pp = -1.780000e+24,
               m_rt =  7.240000e+23,
               m_rp = -2.590000e+23,
               m_tp =  1.390000e+24
               )

# Check what the current Mw is
mt = list(mt_dict.values())
sum_mt_squared = sum([_ ** 2 for _ in mt])
m0 = (1 / sqrt(2)) * (sqrt(sum_mt_squared))
mw = (2/3) * log10(m0) - 10.7
print(f"current mw = {mw:.2f}")
beachball(mt, size=100, linewidth=2, facecolor="b")

# Kludgey manual scale factor
scale_factor = 15.5
mt = [_ * scale_factor for _ in mt]
sum_mt_squared = sum([_ ** 2 for _ in mt])
m0 = (1 / sqrt(2)) * (sqrt(sum_mt_squared))
mw = (2/3) * log10(m0) - 10.7
print(f"scaled mw = {mw:.2f}")
beachball(mt, size=100, linewidth=2, facecolor="r")

for i, (key, val) in enumerate(mt_dict.items()):
    suffix = key.split("_")[1]
    print(f"M{suffix}:      {val * scale_factor:13.6e}")

