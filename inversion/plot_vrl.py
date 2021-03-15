import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

event = None
station = None
comp = None

vrl = np.loadtxt(sys.argv[1], dtype=str, delimiter=", ")
vrl_dict = {"event": vrl[:,0], "station": vrl[:,1], "comp": vrl[:,2], 
            "vrl": vrl[:,3].astype(float)}
df = pd.DataFrame(vrl_dict)
if event:
    df = df.loc[df["event"] == event]
if station:
    df = df.loc[df["station"] == station]
if comp:
    df = df.loc[df["comp"] == comp]

data = df["vrl"].to_numpy()

# Weights makes the y-axis percentage rather than count
plt.hist(data, np.arange(-20, 21, 1), weights=np.ones(len(data)) / len(data),
         edgecolor="k", zorder=5, color="steelblue", linewidth=1.5)
plt.xlabel("Waveform VRL")
plt.ylabel(f"Percentage (N={len(data)})")
plt.title(f"Whole Waveform VRL 6-30s\n"
          f"event={event}; station={station}; comp={comp}\n"
          f"mean={data.mean():.2f}; min={data.min():.2f}; max={data.max():.2f}")
plt.axvline(x=0, linestyle="--", c="k", zorder=6)
plt.axvline(x=data.mean(), linestyle="-", c="orange", zorder=6)
plt.xlim([-15, 15])
plt.ylim([0, 0.2])
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
for axis in ["top", "bottom", "left", "right"]:
    plt.gca().spines[axis].set_linewidth(2)
plt.savefig(sys.argv[1].replace("txt", "png"))
