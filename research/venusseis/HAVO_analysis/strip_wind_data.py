"""
Save wind speed and direction as separate MiniSEED files.
Original files provided by Il-Sang
"""
import os
from glob import glob
from obspy import read

path_out = "/Users/chow/Work/research/venusseis/HAVO_data_analysis/wind_data"
for fid in glob("*/*.miniseed"):
    st = read(fid)
    for channel in ["HWD", "HWS"]:
        st_tmp = st.select(channel=channel)
        comp = channel[-1]
        for tr in st_tmp:
            tr.stats.network = "VH"
            tr.stats.location = ""
            _, _, comp = tr.stats.channel 
            tr.stats.channel = f"LW{comp}"

        start = tr.stats.starttime
        fid_out = f"{tr.get_id()}.{start.year}.{start.julday:0>3}"
        print(fid_out)
        st_tmp.write(os.path.join(path_out, fid_out), format="MSEED")

