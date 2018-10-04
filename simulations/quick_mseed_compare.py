"""
During mesh testing, we often compare waveforms produced on different meshes,
to make sure that the new meshes we are generating are producing the same,
or very similar waveforms as those of our other meshes. This is to ensure
that we aren't introducing any type of numerical artifacts that might sneak
into our simulations and mess with our results.

This script simply reads in the available synthetic data (assumed to already
by in .mseed format, if not can be converted using the script ascii2mseed),
filters to the assumed filter bounds of our tomography problem, and plots them
together with the appropriate color and label.
"""
import os
import glob

from obspy import read
import matplotlib.pyplot as plt


A = "NEWRUN"
B = "OLDRUN"

for fid_A in glob.glob(os.path.join("./", A, "*.mseed")):
    fid_B = os.path.join(os.path.join("./", B), os.path.basename(fid_A))
    if not os.path.exists(fid_B):
        continue
    for fid, color_, label_ in zip([fid_A, fid_B], ['r', 'k'], [A, B]):
        st = read(fid)
        st.filter('bandpass', freqmin=1/30, freqmax=1/6)
        plt.plot(st[0].data, color_, label=label_)

    plt.title("{0} {1} vs {2} at 2 to 30s bandpass".format(
        st[0].stats.station, A, B)
        )
    plt.grid()
    plt.legend()
    plt.show()

