import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal import PPSD
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm

# because I already set some colors
color_dict = {"MWZ":"b","BKZ":"orange","KNZ":"g","RTZ":"r","PUZ":"purple",
            "PXZ":"brown","HAZ":"c","URZ":"k","MXZ":"y","TSZ":"pink",
            "OPRZ":"chartreuse"}

# set path
npz_path = "/seis/prj/fwi/bchow/spectral/ppsd_arrays/"
npz_files = glob.glob(npz_path + "*HHZ*001-365.npz")

# start figure
f = plt.figure(dpi=200)
# loop through filenames
for fid in npz_files:
    sta,cha,year = os.path.basename(fid).split(".")[:3]
    ppsd = PPSD.load_npz(fid)
    mean = ppsd.get_mean()
    plt.plot(mean[0],mean[1],label=sta,color=color_dict[sta])

nlnm_x,nlnm_y = get_nlnm()
nhnm_x,nhnm_y = get_nhnm()
plt.plot(nlnm_x,nlnm_y,'gray',alpha=0.7)
plt.plot(nhnm_x,nhnm_y,'gray',alpha=0.7)

plt.legend(ncol=2)
plt.xlim([0.2,100])
plt.xscale("log")
plt.xlabel("Period (s)")
plt.ylabel("Amplitude [m^2/s^4/Hz][dB]")
plt.title("Mean values of year-long PPSD\'s for GEONET permanent seismometers")
plt.grid()

plt.show()
