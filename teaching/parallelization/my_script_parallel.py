from concurrent.futures import ProcessPoolExecutor
import matplotlib.pyplot as plt
from glob import glob
from obspy import read, Stream

def filter_trace(st, i):
    """
    Filters waveform for a given trace with a
    .01-1Hz bandpass
    :type st: obspy.core.stream.Stream
    :param st: waveform stream to index
    :type i: int
    :param i: index used to select trace
    :rtype: obspy.core.trace.Trace
    :return: filtered waveform trace
    """
    tr = st[i]  # copy data of station
    tr.detrend("demean")
    tr.filter(type="bandpass", freqmin=.01,
              freqmax=1)

# Read in waveform data from disk
st = Stream()
for fid in glob("/home/bchow/Work/comms/faculty_seminar/2009-04-07T201255_SOUTHERN_ALASKA/SAC/*.sac"):
   st += read(fid)

# Simple preprocessing and filtering option
with ProcessPoolExecutor() as executor:
	futures = [executor.submit(filter_trace, st, i)
               for i in range(len(st))]

st.detrend("demean")
st.filter(type="bandpass", freqmin=.01, freqmax=1)

# Plot waveform data in a record section
f, ax = plt.subplots()

for i, tr in enumerate(st):
   tr.data /= tr.data.max()  # normalize between -1 and 1
   plt.plot(tr.times(), tr.data + i, c="k", 
linewidth=.5, linestyle="-")

plt.savefig("test.png")
