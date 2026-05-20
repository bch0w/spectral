"""
Remove instrument response from SmartSolo IGU-BD3C-5 Instruments

ObsPy response removal documentation
https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.remove_response.html#obspy.core.trace.Trace.remove_response

.. rubric::

    $ python remresp_smartsolo.py PATH/TO/DATA/*  # or similar file tag

.. requires::

    ObsPy
"""
import sys
import os
from obspy import read, read_inventory

# =========================== USER PARAMETERS ==================================
# IRIS Nominal Response Library (NRL) path for SmartSolo IGU-BD3C-5
# Note: This points toa  specific instrument page with pre-set parameters 
# regarding pre-filtering and DC offset, you may need to change the input link
# if your instrument parameters are different
RESP_LINK = ("https://service.iris.edu/irisws/nrl/1/combine?instconfig="
             "datalogger_DTCC_SmartSolo-IGU-BD3C-5_PD0_FR100_FPLP_DFDC_OUcount"
             "&format=stationxml")
OUTPUT_PATH = "./resp_removed"
PRE_FILT = [.001, .005, 120, 125]  # pre-filter, if 'None', will not apply
OUTPUT = "VEL"  # DISP=displacement, VEL=velocity, ACC=acceleration
# ==============================================================================

# Begin response removal
inv = read_inventory(RESP_LINK)
if not os.path.exists(OUTPUT_PATH):
    os.makedirs(OUTPUT_PATH)

files = sys.argv[1:]
print(f"removing response for {len(files)} waveforms")

for fid in files:
    fid_base = os.path.basename(fid)
    print(fid_base, end="... ")

    # Check if this data has already been processed
    path_out = os.path.join(OUTPUT_PATH, fid_base)
    if os.path.exists(path_out):
        print("skipped, already processed")
        continue

    # Read data, get metadata
    try:
        st = read(fid)
    except (TypeError, IsADirectoryError):
        print("skipped, unknown file format")
        continue

    net, sta, loc, cha = st[0].get_id().split(".")

    # Rebuild the inventory by name assuming only trace per file
    inv[0].code = net
    inv[0][0].code = sta
    inv[0][0][0].code = cha
    inv[0][0][0].location_code = loc

    # Remove response with optional options
    st.remove_response(inventory=inv, pre_filt=PRE_FILT, output=OUTPUT)

    # Write out new file with response removed
    st.write(path_out, format="MSEED")
    print("done")

