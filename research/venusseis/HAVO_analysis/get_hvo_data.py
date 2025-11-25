"""
Grab HVO data from IRIS
"""
import os
from glob import glob
from obspy import read, UTCDateTime, Stream
from obspy.clients.fdsn import Client


c = Client("IRIS")
path_out = "/Users/chow/Data/nodes/2025-03_HAVO_VENUS/hvodata"

# Gather for each station0
for code in ["HV.DESD..EH?", "HV.MITD.*.?H?", "HV.KAED..EH?", "IU.POHA.*.BH?"]:
    net, sta, loc, cha = code.split(".")

    # Gather for each day of the deployment
    # for julday in range(81, 124, 1):  # temp deploy length
    for julday in range(87, 88, 1):  # temp deploy length
        start = UTCDateTime(f"2025-{julday:0>3}T00:00:00") 
        end = UTCDateTime(f"2025-{julday:0>3}T23:59:59.59999") 

        # Check if any file exists
        fid_check = f"{net}.{sta}.{loc}.{cha}.{start.year}.{start.julday:0>3}"
        if glob(os.path.join(path_out, fid_check)):
            print(f"{fid_check} file exists, skipping")
            continue

        # Raw data
        st = c.get_waveforms(network=net, station=sta, location=loc, 
                             channel=cha, starttime=start, 
                             endtime=end)
        inv = c.get_stations(network=net, station=sta, location=loc, 
                             channel=cha, starttime=start, 
                             endtime=end, level="response")

        # Remove response
        st.remove_response(inventory=inv, output="VEL")

        # Write to disk
        for tr in st:
            fid_tr = f"{tr.id}.{start.year}.{start.julday:0>3}"
            print(fid_tr)
            tr_path = os.path.join(path_out, fid_tr)
            tr.write(tr_path, "MSEED")
        del st


# for comp in ["E", "N", "Z"]:
#     # Check existing files
#     fid_check = f"{net}.{sta}.*.??{comp}.{start.year}.{start.julday}"
#     check_bool = glob(os.path.join(path_out, fid_check))
#     if bool(check_bool):
#         print(f"{fid_check} exists, skipping")
#         continue
