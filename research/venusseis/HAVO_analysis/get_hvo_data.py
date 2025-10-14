"""
Grab HVO data from IRIS
"""
import os
from glob import glob
from obspy import read, UTCDateTime, Stream
from obspy.clients.fdsn import Client


c = Client("IRIS")
net = "HV"
path_out = "/Users/chow/Work/research/venusseis/HAVO_data_analysis/node_data/hvodata"

# Gather for each station
for sta in ["DESD", "MITD", "KAED"]:
    # Gather for each day of the deployment
    for julday in range(81, 124, 1):  # temp deploy length
        start = UTCDateTime(f"2025-{julday:0>3}T00:00:00") 
        end = UTCDateTime(f"2025-{julday:0>3}T23:59:59.59999") 
        # end = UTCDateTime(f"2025-{julday:0>3}T01:00:00.")   # TEST

        # Raw data
        for comp in ["E", "N", "Z"]:
            # Check existing files
            fid_check = f"{net}.{sta}.*.??{comp}.{start.year}.{start.julday}"
            check_bool = glob(os.path.join(path_out, fid_check))
            if bool(check_bool):
                print(f"{fid_check} exists, skipping")
                continue

            st = c.get_waveforms(network=net, station=sta, location="*", 
                                 channel=f"EH{comp}", starttime=start, 
                                 endtime=end)

            for tr in st:
                fid_tr = f"{tr.id}.{start.year}.{start.julday:0>3}"
                print(fid_tr)
                tr_path = os.path.join(path_out, fid_tr)
                breakpoint()
                st.write(tr_path, "MSEED")
            del st


