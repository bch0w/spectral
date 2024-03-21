"""
Get station pair information from the CC stacks based on file names
"""
import os 
from glob import glob
from obspy import read


netstas_a = []
netstas_b = []
lats = dict()
lons = dict()
with open("station_pairs.txt", "w") as f:
    for dir_ in glob("SAC_I3_stack_4_Zendo/*_stack"):
        phase, _, typ, _ = os.path.basename(dir_).split("_")
        for dir_ in sorted(glob(os.path.join(dir_, "*_*"))):
            for i, fullpath in enumerate(sorted(glob(os.path.join(dir_, "*")))):
                fid = os.path.basename(fullpath)
                fid = os.path.splitext(fid)[0]
                _, net_a, sta_a, net_b, sta_b = fid.strip().split("_")
                if i == 0:
                    f.write(f"{phase.upper()}  {typ.upper()}  {net_a}.{sta_a}\n")
                f.write(f"\t{net_b}.{sta_b}\n")

                netsta_a = f"{net_a} {sta_a}"
                netsta_b = f"{net_b} {sta_b}"
                if netsta_a not in netstas_a:
                    netstas_a.append(netsta_a)
                if netsta_b not in netstas_b:
                    netstas_b.append(netsta_b)
                    st = read(fullpath)
                    lats[netsta_b] = st[0].stats.sac.stla
                    lons[netsta_b] = st[0].stats.sac.stlo


assert(len(netstas_a) == len(netstas_b))

with open("STATIONS_LIU2022", "w") as f:
    for netsta in sorted(netstas_b):
        net, sta = netsta.split(" ")
        f.write(f"{sta:>6}{net:>6}{lats[netsta]:11.4f}{lons[netsta]:11.4f}{0:7.1f}{0:7.1f}\n")

