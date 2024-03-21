"""
Generate a waveform that looks like an 'm' if you squint
"""
import numpy as np
from scipy.signal import tukey as window
import matplotlib.pyplot as plt

tail = 2.75
x = np.linspace(-1 * tail * np.pi, tail * np.pi, 1000)
# Plot a shifted version for that stereo effect
x2 = x - np.pi*.2
y = -1 * np.cos(x)
taper = window(len(x), alpha=.8)
y_win = np.multiply(y, taper)

c_syn = "orangered"
lw=22
plt.plot(x[:970], y_win[:970], c=c_syn, linewidth=lw-5)
plt.plot(x2, y_win, c="k", linewidth=lw)
plt.axis("off")
plt.savefig("m_in_atom.png", transparent=True)
