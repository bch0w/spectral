import os
from pyatoa import Manager
from pyasdf import ASDFDataSet as asdf
import matplotlib.pyplot as plt

periods = [3, 4]

with asdf("2015p768477_birch.h5") as ds:
    for station in ds.waveforms.list():
        f, axes = plt.subplots(3)
        mgmt = Manager(ds=ds)
        try:
            for i, min_period in enumerate(periods):
                mgmt.load(station, path="i11/s05")
                mgmt.config.min_period = min_period
                mgmt.standardize().preprocess()
                for j, tr in enumerate(mgmt.st_obs):
                    axes[j].plot(tr.times(), tr.data, 
                                 c={0: "k", 1: "r"}[i])
                mgmt.reset()

            axes[0].set_title(station)
            axes[1].set_ylabel("Displacement [m]")
            axes[2].set_xlabel("Time [s]")
            fid = f"{station}_2015p768477_{periods[0]}vs{periods[1]}s.png"

            plt.savefig(os.path.join("figures", fid), dpi=100)
            plt.close()
        except Exception as e:
            print(f"{station}: {e}")
            plt.close()
            continue

