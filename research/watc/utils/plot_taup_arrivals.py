"""
Use TauPy to print out arrival times for a given source depth and station 
distance in degrees
"""
import sys
import matplotlib.pyplot as plt
import numpy as np
from pysep import read_sem
from obspy.taup import TauPyModel, plot_travel_times
from obspy.geodetics import kilometers2degrees


fid = sys.argv[1]
st = read_sem(fid, source="CMTSOLUTION", stations="STATIONS")
distance_in_degrees = kilometers2degrees(kilometer=st[0].stats.sac.dist)
source_depth_in_km = st[0].stats.sac.evdp

# Model parameters
model = TauPyModel(model="prem")
# phase_list = ["ttbasic"]
# phase_list = ["P", "PP", "pS", "PS", "Pdiff", "pPdiff", "PKP", "PcP",
#               "S", "SS", "SP", "sS", "Sdiff", "sSdiff", "SKS", "ScP"]
phase_list = ["P", "S", "PP", "SS"]

arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km,
                                  distance_in_degree=distance_in_degrees,
                                  phase_list=phase_list
                                  )

# Preprocess waveform
st.taper(max_percentage=0.05)
st.differentiate()
st.filter("lowpass", freq=1/30)

# Plot the waveform
f, ax = plt.subplots(figsize=(8, 4), dpi=200)

# Custom get time axis
t_max = 2300
x_max = int(t_max // st[0].stats.delta)
t = st[0].times()[:x_max]
d = st[0].data[:x_max]

if False:
    t = st[0].times()
    d = st[0].data

plt.plot(t, d, color="k", lw=1)
plt.xlabel("Time [s]")
# plt.xlim([0, st[0].times().max()])
plt.xlim([0, t_max])
plt.ylabel("Velocity [m/s]")

# Set major/minor tick marks
major_ticks = np.arange(t.min(), t.max(), 250)
minor_ticks = np.arange(t.min(), t.max(), 50)
ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor=True)

# Plot phase arrivals
n = len(arrivals)
maxamp = d.max()
yvals = np.linspace(maxamp/4, maxamp, n)
for i, arrival in enumerate(arrivals):
    time_s = arrival.time
    name = arrival.name
    plt.axvline(x=time_s, c=f"C{i}")
    plt.text(x=time_s, y=yvals[i], s=name)

plt.savefig("waveform.png")

plt.close("all")

# Plot arrivals
arrivals = model.get_ray_paths(source_depth_in_km=source_depth_in_km,
                               distance_in_degree=distance_in_degrees,
                               phase_list=phase_list,)

arrivals.plot_rays(phase_list=phase_list, legend=True, show=False)
plt.savefig("arrivals.png")
plt.close("all")

# Plot Traveltimes
plot_travel_times(source_depth=source_depth_in_km, phase_list=phase_list, 
                  show=False)
plt.savefig("travel_times.png")
