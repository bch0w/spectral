"""
Testing misfit quantification of EGC from Liu et al. (2022) and SGF generated
by SPECFEM3D_GLOBE (w/ S20RTS and Crust2.0). Generated for master station
IU.COLA using three individual Z, N and E point force simulations. N and E 
components were rotated to TT and RR (not used) using the rotation matrix 
from Wang et al. (2019)

Comparing ZZ and TT components, which are Rayleigh and Love, respectively
"""
import os
import matplotlib.pyplot as plt

from glob import glob
from obspy import read
from pysep.utils.io import read_sem

min_period = 10
max_period = 20
master_station="IU_COLA"

# Loop through each phase (rayleigh, love) and type (hyp, ell)
for stack_path in sorted(glob("EGF/*_stack")):
    stack_fid = os.path.basename(stack_path)
    phase, _, type_, *_ = stack_fid.split("_")

    # Determine the components based on the seismic wave used
    if phase == "Ray":
        comp = "Z"
    elif phase == "Lov":
        comp = "T"

    # Loop through available EGF
    for fullpath in sorted(glob(os.path.join(stack_path, master_station, "*"))):
        fid = os.path.basename(fullpath)
        # ex. I3_IU_COLA_TA_L19K.SAC
        _, _, _, net, sta = fid.split(".")[0].split("_")

        # Find matching SGF
        sgf_path = os.path.join(f"{comp}{comp}_SGF", 
                                f"{sta}.{net}.BX{comp}.sem.ascii")
        # Not all SGF exist based on the smaller synthetic domain
        if not os.path.exists(sgf_path):
            continue

        st_egf = read(fullpath) 

        # Change the stats so the figure looks cleaner
        st_egf[0].stats.station = sta
        st_egf[0].stats.location = "EGF"
        st_egf[0].stats.channel = f"{comp}{comp}"

        st_sgf = read_sem(sgf_path, origintime=st_egf[0].stats.starttime)  
        st_sgf[0].stats.location = "SGF"

        st = st_egf + st_sgf

        # Data starts first but synthetics end first
        st.trim(st[0].stats.starttime, st[1].stats.endtime)

        st.taper(max_percentage=0.05)
        st.filter("bandpass", freqmin=1/max_period, freqmax=1/min_period)

        # Normalize to 1
        for tr in st:
            tr.data /= tr.data.max()

        f, axs = plt.subplots(3, sharex=True)
        f.suptitle(f"{master_station.replace('_', '.')} -> {net}.{sta} "
                   f"{comp}{comp} {min_period}-{max_period}s [{phase} {type_}]")

        # Plot the EGF by itself
        axs[0].plot(st[0].times(), st[0].data, "k", linewidth=1., label="EGF")

        # Plot the SGF by itself
        axs[1].plot(st[1].times(), st[1].data, "r", linewidth=1., label="SGF")

        # Plot the two together
        axs[2].plot(st[0].times(), st[0].data, "k", linewidth=1., label="EGF")
        axs[2].plot(st[1].times(), st[1].data, "r", linewidth=1., label="SGF")


        axs[2].set_xlabel("Time [s]")
        axs[1].set_ylabel("Normalized Amplitude")

        plt.savefig(f"figures_10-20s/{net}_{sta}_{phase.lower()}_{type_}.png")
        plt.close()
                

        






