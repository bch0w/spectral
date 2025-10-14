import os
from glob import glob
from obspy import read



path_ = "/Users/chow/Data/nodes/2025-03_HAVO_VENUS/"
for sta_path in glob(os.path.join(path_, "*_*")):
    sta_name_full = os.path.basename(sta_path)  #  e.g., DESD_bare
    sta_name, sta_type = sta_name_full.split("_")
    sta_loc = {"bare": "EX", "bucket": "CV", "windshield": "WS"}[sta_type]
    for full_path in glob(os.path.join(sta_path, "*.miniseed")):
        _fid = os.path.basename(full_path)

        st = read(full_path)
        st[0].stats.network = "VH"  # Venus Hawaii
        st[0].stats.station = sta_name.upper()
        st[0].stats.location = sta_loc
        
        origintime = st[0].stats.starttime
        year = origintime.year
        jday = origintime.julday

        fid_out = f"{st[0].get_id()}.{year}.{jday:0>3}"
        full_path_out = os.path.join(path_, fid_out)

        print(f"{_fid} -> {fid_out}")
        st.write(full_path_out, format="MSEED")
        # os.rename(full_path, full_path_out)

        

