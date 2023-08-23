"""
Creates a directory of symlinks which point to the Liu et al. 2022 EGF data
in a format that SeisFlows can determine
"""
import os
from glob import glob
from pysep.utils.io import read_stations


# input_dir = "/import/c1/ERTHQUAK/bhchow/data/egfs/SAC_I3_stack_4_Zendo"
# output_dir = "/import/c1/ERTHQUAK/bhchow/data/egfs/NALASKA_EGF"
# stations_file = "/import/c1/ERTHQUAK/bhchow/work/akatom/northern_alaska/SFDATA/STATIONS"
input_dir = "/home/bchow/Work/data/egfs/SAC_I3_stack_4_Zendo/"
output_dir = "/home/bchow/Work/data/egfs/NALASKA_EGF"
stations_file = "/home/bchow/Work/work/seisflows/anat_trial/workdir/DATA/STATIONS"

# Used for cutting out stations that are not in the domain
inv = read_stations(stations_file)

for stack in glob(os.path.join(input_dir, "*")):
    phase, _, type_, _ = os.path.basename(stack).split("_")

    kernel = {"Lov": "TT", "Ray": "ZZ"}[phase]
    comp = kernel[-1]

    for src in glob(os.path.join(stack, "*")):
        # Dictate that source station need to be in the simulation domain
        net_src, sta_src = os.path.basename(src).split("_")
        if not bool(inv.select(network=net_src, station=sta_src)):
            continue
        for fid in glob(os.path.join(src, "*")):
            # Separate out the station names from rest of file name
            _fid = os.path.splitext(os.path.basename(fid))[0]
            _, net_src, sta_src, net_rcv, sta_rcv = _fid.split("_")

            # Dictate that rcv station need to be in the simulation domain
            if not bool(inv.select(network=net_rcv, station=sta_rcv)):
                continue

            # Determine the correct save location
            path_out = os.path.join(output_dir, type_.lower(), 
                                    f"{net_src}_{sta_src}", kernel)
            if not os.path.exists(path_out):
                os.makedirs(path_out)
        
            # Label the station correctly, L because it's 1Hz data, X because
            # it's a time series derived from observational data and then the
            # corresponding component
            filename_out = f"{net_rcv}.{sta_rcv}.LX{comp}.SAC"
            fid_out = os.path.join(path_out, filename_out)
            
            if not os.path.exists(fid_out):
                os.symlink(fid, fid_out)
