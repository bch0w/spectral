"""
Convert USGS Stream Gauge data, whichi is provided in ASCII format, into a 
an ObsPy Stream and save as a MiniSEED file so that it can be plotted easily
"""
import sys
import os
import numpy as np


# Input Parameters
path = sys.argv[1]
units = "m"  # ft, in, m, cm
relative = "min"  # 'origin', 'min', 'mean', 'med', 'abs'


# Begin processing
assert(os.path.exists(path))

# Get data
data = np.loadtxt(path, skiprows=28, usecols=[2,4], delimiter="\t", dtype=str)
time, amplitude = data.T
amplt

# Convert units
if units == "ft":
    height = height_ft
elif units == "in":
    height = height_ft * 12
elif units == "m":
    height = np.array([_ * 0.3048 for _ in height_ft.astype(float)])
elif units == "cm":
    height = np.array([_ * 30.48 for _ in height_ft.astype(float)])
else:
    print("`units` should be in: ft, in, m, cm")
    sys.exit()

# Subset data where we are plotting waveforms to get the correct ylims
idx = np.where((times > xvals.min()) & (times < xvals.max()))

if args.sg_rel:
    height -= height[idx].min()

# Plot on the same axis as the waveform
twax = ax.twinx()

# # Plot the gradient, not very helpful
# dx = (times[1] - times[0]) * 86400
# grad_height = np.gradient(height[idx], dx)
# twax.plot(times[idx], grad_height, "o-", lw=1, c="C0", 
#         label="Phelan Cr. Gage", zorder=5, markersize=1.25, 
#         alpha=0.5)

twax.plot(times[idx], height[idx], "o-", lw=1, c="C0", 
            label="Phelan Cr. Gage", zorder=5, markersize=1.25, 
            alpha=0.5)

_ylabel = f"Stream Height [{args.sg_units}]"
if args.sg_rel:
    _ylabel = f"Relative {_ylabel}"
twax.set_ylabel(_ylabel, rotation=-90, labelpad=20)
