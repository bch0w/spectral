"""
Small utility function to change network codes in MSEED data from XX to 2P
for the BEACON deployment.
"""
import os
from glob import glob
from obspy import read


path = "/home/bchow/Work/data/mseed/BEACON"

for year in glob(os.path.join(path, "201?")):
    for net in glob(os.path.join(year, "??")):
        for sta in glob(os.path.join(net, "RD??")):
            for cha in glob(os.path.join(sta, "HH?.D")):
                for fid in glob(os.path.join(cha, "XX.*")):
                    fid_new = fid.replace("XX", "2P")
                    print(f"{os.path.basename(fid)} -> {os.path.basename(fid_new)}")
                    st = read(fid)
                    for tr in st:
                        tr.stats.network = "2P"
                    st.write(fid_new, format="MSEED")
                    os.remove(fid)


