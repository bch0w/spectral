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

from pyatoa import Config, Manager
from pyatoa.core.manager import ManagerError

i=0
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
        st_egf[0].data /= st_egf[0].data.max()

        # Strip the SAC header because it causes Pyflex to fail
        del st_egf[0].stats.sac

        st_sgf = read_sem(sgf_path, origintime=st_egf[0].stats.starttime)  
        st_sgf[0].stats.location = "SGF"
        st_sgf[0].data /= st_sgf[0].data.max()

        # Start Pyatoa workflow
        cfg = Config(min_period=8, max_period=80, st_obs_type="syn", 
                     st_syn_type="syn", component_list=[comp],
                     adj_src_type="cc_traveltime")

        mgmt = Manager(config=cfg, st_obs=st_egf, st_syn=st_sgf)
        try:
            mgmt.standardize()
            mgmt.preprocess()
            mgmt.window()
            mgmt.measure()
            mgmt.plot(show=False, 
                      save=f"pyatoa_figures/{net}_{sta}_{phase.lower()}_{type_}.png")
        except ManagerError:
            continue

                

        






