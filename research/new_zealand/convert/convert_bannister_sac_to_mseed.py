"""
Stephen Bannister provided his temporary waveform data in SAC format separated 
by event ID. Pyatoa needs this data in SEED format so convert to mseed and
rearrange in the correct directory structure
"""
import os
from glob import glob
from obspy import read

# Data dir and output
data_dir = "/home/chowbr/gns03247/project/yoshi/BANNISTER/{dep}/*/*"
out_dir = "/home/chowbr/gns03247/project/bchow/data/waveforms"

for dep in ["HD", "GA"]:
    for fid in glob(data_dir.format(dep=dep)):
        st = read(fid)
        year = str(st[0].stats.starttime.year)
        julday = str(st[0].stats.starttime.julday)
        code = st[0].get_id()
        net, sta, loc, cha = code.split(".")

        dir_structure = os.path.join(out_dir, dep, year, net, sta, f"{cha}.D")
        if not os.path.exists(dir_structure):
            os.makedirs(dir_structure)

        fid_out = f"{net}.{sta}.{loc}.{cha}.D.{year}.{julday}"

        print(f"{os.path.basename(fid)} -> {fid_out}")

        st.write(os.path.join(dir_structure, fid_out), format="MSEED")



