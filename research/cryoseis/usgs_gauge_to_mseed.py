"""
Convert USGS Stream Gauge data, whichi is provided in ASCII format, into a 
an ObsPy Stream and save as a MiniSEED file so that it can be plotted easily
"""
import sys
import os
import numpy as np
from obspy import UTCDateTime, Stream, Trace

# Input Parameters
path = sys.argv[1]
units = "cm"  # ft, in, m, cm
relative = 1460.83  # 'origin', 'min', 'mean', 'med', 'abs', None or float
time_offset = +8  # local time to UTC [hours]

# Begin processing
assert(os.path.exists(path))

# Get data
data = np.loadtxt(path, skiprows=28, usecols=[2, 6], delimiter="\t", dtype=str)
time, amplitude = data.T  # amplitude is discharge in gauge hheight in feet

# Create time axis
starttime = UTCDateTime(time[0]) 
endtime = UTCDateTime(time[-1])

dt = UTCDateTime(time[1]) - starttime
tlen = endtime - starttime

starttime += (time_offset * 60 * 60)  # time offset hours -> seconds

# Convert units from feet to user-selected
height = amplitude.astype(float)

# Convert units 
if units == "ft":
    pass
elif units == "in":
    height *= 12
elif units == "m":
    height *= 0.3048
elif units == "cm":
    height *= 30.48
else:
    print("`units` should be in: ft, in, m, cm")
    sys.exit()

# Set values to relative
if relative == "min":
    height -= height.min()
elif relative == "origin":
    height -= height[0]
elif relative == "mean":
    height -= height.mean()
elif relative == "abs":
    height = np.abs(height)
elif isinstance(relative, (float,int)):
    height -= relative
else:
    pass


# Create an ObsPy Stream
tr = Trace(data=height)
tr.stats.delta = dt
tr.stats.starttime = starttime
tr.stats.network = "US"
tr.stats.station = "PHCR"
tr.stats.location = ""
tr.stats.channel = "USZ"

st = Stream(traces=[tr])
st.write("US.PHCR..USZ", format="MSEED")


