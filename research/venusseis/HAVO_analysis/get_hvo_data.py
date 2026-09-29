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
# codes = ["HV.DESD..EH?", "HV.MITD.*.?H?", "HV.KAED..EH?", "IU.POHA.*.BH?"]
codes = ["HV.DESD..EHZ"]
juldays = [87]  # range(81, 124, 1)

# Gather for each station0
for code in codes:
    net, sta, loc, cha = code.split(".")

    # Gather for each day of the deployment
    for julday in juldays:
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
        st.remove_response(inventory=inv, output="VEL", pre_filt=pre_filt)

        # Write to disk
        for tr in st:
            fid_tr = f"{tr.id}.{start.year}.{start.julday:0>3}"
            print(fid_tr)
            tr_path = os.path.join(path_out, fid_tr)
            tr.write(tr_path, "MSEED")
        del st


