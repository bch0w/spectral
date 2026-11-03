"""
Grab HVO data from IRIS
"""
import os
from glob import glob
from obspy import read, UTCDateTime, Stream
from obspy.clients.fdsn import Client


c = Client("EARTHSCOPE")
path_out = "./"
pre_filt = [.001, .005, 120, 125]
codes = ["HV.DESD..EH?", "HV.MITD.*.?H?", "HV.KAED..EH?", "IU.POHA.*.BH?"]
starttime = UTCDateTime("2025-087T06:00:00")
endtime = starttime + 60 * 60 * 3

# Gather for each station0
for code in codes:
    net, sta, loc, cha = code.split(".")

    # Gather for each day of the deployment
    # Check if any file exists
    fid_check = f"{net}.{sta}.{loc}.{cha}.{starttime.year}.{starttime.julday:0>3}"
    if glob(os.path.join(path_out, fid_check)):
        print(f"{fid_check} file exists, skipping")
        continue

    # Raw data
    st = c.get_waveforms(network=net, station=sta, location=loc, 
                         channel=cha, starttime=starttime, 
                         endtime=endtime)
    inv = c.get_stations(network=net, station=sta, location=loc, 
                         channel=cha, starttime=starttime, 
                         endtime=endtime, level="response")

    # Remove response
    st.remove_response(inventory=inv, output="VEL", pre_filt=pre_filt)

    # Write to disk
    for tr in st:
        fid_tr = f"{tr.id}.{starttime.year}.{starttime.julday:0>3}"
        print(fid_tr)
        tr_path = os.path.join(path_out, fid_tr)
        tr.write(tr_path, "MSEED")
    del st


