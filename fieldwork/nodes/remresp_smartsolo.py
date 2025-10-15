"""
Brute force remove response from SmartSolo Instruments

- Preamp (dB): 0 dB(factory default, rarely changed)
- Final Sample Rate: 100 Hz (10 ms)
- Final Filter Phase: Linear phase(LP)
- DC Filter: Off (Low Cut Filter: Disabled)
"""
import sys
import os
from obspy import read, read_inventory

# IRIS NRL for 
inv = read_inventory("https://service.iris.edu/irisws/nrl/1/combine?instconfig="
                     "datalogger_DTCC_SmartSolo-IGU-BD3C-5_"
                     "PD0_FR100_FPLP_DFDC_OUcount&format=stationxml")

path_out = "/Users/chow/Work/research/venusseis/HAVO_data_analysis/node_data/remresp"

for fid in sys.argv[1:]:
    fid_out = os.path.basename(fid)
    full_out = os.path.join(path_out, fid_out)
    if os.path.exists(full_out):
        continue

    st = read(fid)
    net, sta, loc, cha = st[0].get_id().split(".")

    # Rebuild the inventory by name assuming only one of each
    inv[0].code = net
    inv[0][0].code = sta
    inv[0][0][0].code = cha
    inv[0][0][0].location_code = loc

    st.remove_response(inventory=inv)
    print(fid_out)
    st.write(full_out, format="MSEED")

