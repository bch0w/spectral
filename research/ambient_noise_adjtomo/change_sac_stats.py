"""
Reset the origintime and set trace stats for all the ambient noise data
"""
import os
from glob import glob
from obspy import read, UTCDateTime

# Set desired starttime here
starttime = UTCDateTime("2000-01-01T00:00:00")


path = "/home/bchow/Work/data/egfs/NALASKA_EGF"
for dir_ in glob(os.path.join(path, "???")):  # ell, hyp
    for sta in glob(os.path.join(dir_, "??_*")):  # NN_SSS(S)
        for kernel in glob(os.path.join(sta, "??")):  # TT or ZZ
            for fid in glob(os.path.join(kernel, "*.SAC")):
                print(f"{fid}")
                filename = os.path.basename(fid)
                try:
                    net, sta, cha, sac = filename.split(".")
                    st = read(fid)
                    if len(st) != 1:
                        print(f"!!!! {filename}")
                        continue
                    for tr in st:
                        tr.stats.network = net
                        tr.stats.station = sta
                        tr.stats.channel = cha
                        tr.stats.starttime = starttime
                        tr.stats.sac.nzyear = 2000
                        tr.stats.sac.evdp = 0.
                    st.write(fid, format="SAC")
                except Exception as e:
                    print(f"!!!! {filename}: {e}")
                    continue

