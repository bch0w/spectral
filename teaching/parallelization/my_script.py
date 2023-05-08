import matplotlib.pyplot as plt
from glob import glob
from obspy import read, Stream

# Read in waveform data from disk
st = Stream()
for fid in glob("*.sac"):
   st += read(fid)

breakpoint()

# Simple preprocessing and filtering option
st.detrend("demean")
st.filter(type="bandpass", freqmin=.01, freqmax=1)

# Plot waveform data in a record section
f, ax = plt.subplots()

for i, tr in enumerate(st):
   tr.data /= tr.data.max()  # normalize between -1 and 1
   plt.plot(tr.times(), tr.data + i, c="k", 
linewidth=.5, linestyle="-")

plt.savefig("test.png")
