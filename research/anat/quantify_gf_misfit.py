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

from pyatoa import Config, Manager, logger
from pyatoa.core.manager import ManagerError

logger.setLevel("DEBUG")

i=0
master_station="IU_COLA"
path_egf = "../EGF"
path_sgf = "../SGF"

# Loop through each phase (rayleigh, love) and type (hyp, ell)
for stack_path in sorted(glob(f"{path_egf}/*_stack")):
    stack_fid = os.path.basename(stack_path)
    phase, _, type_, *_ = stack_fid.split("_")

    # Determine the components based on the seismic wave used
    if phase == "Ray":
        comp = "Z"
        comp_list = ["Z", "N", "E"]
        continue
    elif phase == "Lov":
        comp = "T"
        comp_list = ["T", "Z"]

    # Loop through available EGF
    for fullpath in sorted(glob(os.path.join(stack_path, master_station, "*"))):
        fid = os.path.basename(fullpath)
        # ex. I3_IU_COLA_TA_L19K.SAC
        _, _, _, net, sta = fid.split(".")[0].split("_")

        # Find matching SGF
        sgf_path = os.path.join(path_sgf, f"{comp}", 
                                f"{net}.{sta}.BX{comp}.sem.ascii")
        # Not all SGF exist based on the smaller synthetic domain
        if not os.path.exists(sgf_path):
            continue

        st_egf = read(fullpath) 

        # Change the stats so the figure looks cleaner
        st_egf[0].stats.station = sta
        st_egf[0].stats.location = "EGF"
        st_egf[0].stats.channel = f"{comp}{comp}"

        # Strip the SAC header because it causes Pyflex to fail
        del st_egf[0].stats.sac

        st_sgf = read_sem(sgf_path, origintime=st_egf[0].stats.starttime)  
        st_sgf[0].stats.location = "SGF"

        # Start Pyatoa workflow
        cfg = Config(min_period=30, max_period=50, st_obs_type="syn", 
                     st_syn_type="syn", component_list=[comp],
                     adj_src_type="multitaper", stalta_waterlevel=0.12,
                     tshift_acceptance_level=25.)

        mgmt = Manager(config=cfg, st_obs=st_egf, st_syn=st_sgf)

        try:
            mgmt.standardize()
            mgmt.preprocess()
            mgmt.stats.time_offset_sec = \
                    st_sgf[0].stats.starttime - st_egf[0].stats.starttime

            # Normalize amplitudes because amplitude information not required
            for tr in mgmt.st_obs:
                tr.data /= tr.data.max()
            for tr in mgmt.st_syn:
                tr.data /= tr.data.max()

            mgmt.window()
            mgmt.measure()
            mgmt.plot(show=False, 
                      save=f"figures/{net}_{sta}_{phase.lower()}_{type_}.png")
            # Workaround to get other components to write out as 0s
            mgmt.config.component_list = comp_list
            mgmt.write_adjsrcs(path="./", write_blanks=True)
            # Reset for next processing step 
            mgmt.config.component_list = [comp]
        except ManagerError as e:
            print(e)
            continue

                

        






