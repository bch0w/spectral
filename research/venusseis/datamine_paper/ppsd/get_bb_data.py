"""
Get raw HVO data from EarthScope
"""
import os
import sys
from glob import glob
from obspy import read, UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNNoDataException


c = Client("EARTHSCOPE")
path_out = "./DATA"
choice = 2

# HVO
if choice == 0:
    if sta in ["DESD", "KAED"]:
        band = "E"
    elif sta in ["MITD"]:
        band = "H"
    else:
        raise NotImplementedError
    year = 2025
    start_day = 85
    end_day = 121
    # codes = ["HV.DESD..EH?", "HV.MITD.*.?H?", "HV.KAED..EH?", "IU.POHA.*.BH?"]
    codes = [f"HV.{sta}.*.{band}H?"]

# AEC Phase 1
elif choice == 1:
    year = 2025
    start_day = 190
    end_day = 206
    codes = [f"AK.MDM.*.BH?", f"AK.H23K.*.BH?"]

# AEC Phase 2
elif choice == 2:
    year = 2025
    start_day = 245
    end_day = 277
    codes = [f"AK.PAX.*.BH?", f"AK.K24K.*.BH?", f"AK.RND.*.BH?"]

if not os.path.exists(path_out):
    os.makedirs(path_out)

# Gather for each station0
for code in codes:
    net, sta, loc, cha = code.split(".")

    # Raw data
    for jday in range(start_day, end_day):
        starttime = UTCDateTime(f"{year}-{jday:0>3}T00:00:00.0000")
        endtime = UTCDateTime(f"{year}-{jday:0>3}T23:59:59.9999")

        # Check if any file exists
        fid_check = f"{code}.{starttime.year}.{starttime.julday:0>3}"
        if glob(os.path.join(path_out, fid_check)):
            print(f"{fid_check} file exists, skipping")
            continue
        try:
            st = c.get_waveforms(network=net, station=sta, location=loc, 
                                 channel=cha, starttime=starttime, 
                                 endtime=endtime)
        except FDSNNoDataException as e:
            print(f"no data for: {fid_check}")
            continue

        components = set([tr.stats.component for tr in st])

        # Write each component separately. Write as a stream incase streams
        # are gappy
        for comp in components:
            st_write = st.select(component=comp)
            fid_st = f"{st_write[0].id}.{starttime.year}.{starttime.julday:0>3}"
            print(fid_st)
            fid_out = os.path.join(path_out, fid_st)
            st_write.write(fid_out, "MSEED")

    # Get a single response file for that whole time period
    starttime = UTCDateTime(f"2025-{start_day:0>3}T00:00:00")
    endtimetime = UTCDateTime(f"2025-{end_day:0>3}T00:00:00")
    inv = c.get_stations(network=net, station=sta, location=loc, 
                         channel=cha, starttime=starttime, 
                         endtime=endtime, level="response")

    inv.write(f"XML/{net}.{sta}.xml", format="STATIONXML")




