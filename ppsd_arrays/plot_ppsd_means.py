import numpy as np
import matplotlib.pyplot as plt
from obspy.signal import PPSD

ppsd_1 = PPSD.load_npz('./BKZ.HHZ.2015.001-365.npz')
ppsd_2 = PPSD.load_npz('./KNZ.HHZ.2015.001-365.npz')

m1 = ppsd_1.get_mean()
m2 = ppsd_2.get_mean()

f = plt.figure()
plt.plot(m1[0],m1[1],'r',label='BKZ (inland)')
plt.plot(m2[0],m2[1],'k',label='KNZ (coast)')
plt.legend()

plt.xlim([0.1,100])
plt.xscale("log")
plt.xlabel("Period (s)")
plt.ylabel("Amplitude [m^2/s^4/Hz][dB]")
plt.title('Mean values of year-long PPSD (2015) for KNZ (coast) and BKZ (inland)')
plt.grid()

plt.show()
