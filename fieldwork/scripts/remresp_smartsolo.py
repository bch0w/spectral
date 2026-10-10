"""
!!!
NOTE PLEASE READ 10/8/26: 
I have combined all functionality of this script and my fairfield
response removal script into `remresp.py`. Please use that script moving forward
I will no longer update or maintain this script but I will leave it in place
for reference
!!!

Remove instrument response from SmartSolo IGU-BD3C-5 Instruments using
EarthScope Nominal Response Library URLs (NRL)

ObsPy response removal documentation
https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.remove_response.html#obspy.core.trace.Trace.remove_response

.. rubric::
    
    Select user-parameters based on datalogger configuration below and then
    run from the command line as follows:
    
    $ python remresp_smartsolo.py PATH/TO/DATA/XX.*  # or similar file tag

.. requires::

    ObsPy

.. kludge::

    The warning message below is given if our input data is in units 'mV' 
    because this does not match ObsPy's internal mapping dictionaries. This is
    fine because we only want the gain that is added from this stage so we 
    rename the intput units. If the output units are counts then this should
    not be an issue

    UserWarning: The unit 'MV' is not known to ObsPy. It will be passed in to 
    evalresp as 'undefined'. This should result in evalresp using the response 
    as is, without adding any integration or differentiation and the 'output'
    parameter (here: 'VEL') not having any effect.`
"""
import sys
import os
from obspy import read, read_inventory

# =========================== USER PARAMETERS ==================================
# Response building parameters. See link below for options (case-sensitive):
# https://ds.iris.edu/ds/nrl/datalogger/dtcc/smartsolo-igu-bd3c-5/
PREAMP_DB = 0  # 0 or 6
SAMPLE_RATE = 100  
FILTER_PHASE = "LP"  # LP (linear phase), MP (minimum phase)
DC_FILTER = "Off"  #  1, DC, Off
OUTPUT_UNITS = "mV"  # count, mV

# Response removal parameters
PRE_FILT = [.001, .005, 120, 125]  # pre-filter, if 'None', will not apply
OUTPUT = "VEL"  # output unit; DISP=displacement, VEL=velocity, ACC=acceleration

#  Output parameters
OUTPUT_PATH = "./resprmv"  # path to save new data files with resp. removed
OVERWRITE = True  # if False, will not overwrite existing files

# Rename internal stream stats to match filename. Filename MUST be in the 
# following format or this will fail: NN.SSSS.LL.CCC.YYYY.JJJ
# where N=network, S=station, L=location, C=channel, Y=year, J=julian day
RENAME = True 
# ==============================================================================

# Build EarthScope NRL link based on the user-set parameters
# We assume these are the BD3C-5 broadband instruments
RESP_LINK = (
        "https://service.earthscope.org/irisws/nrl/1/combine?"
        "instconfig=sensor_DTCC_DT-SOLO-BB_LP5_SG209.4_STgroundVel:"
        "datalogger_DTCC_SmartSolo-IGU-BD3C-5_"
        f"PD{PREAMP_DB}_FR{SAMPLE_RATE}_FP{FILTER_PHASE}_DF{DC_FILTER}_OU{OUTPUT_UNITS}&"
        "format=stationxml&nodata=404"
        )
print(f"response url: {RESP_LINK}")
try:
    inv = read_inventory(RESP_LINK)
except Exception as e:
    print(e)
    sys.exit(-1)
print(inv[0][0][0].response)

# See `Kludge` comment in top docstring for explanation of this operation
if OUTPUT_UNITS == "mV":
    resp = inv[0][0][0].response
    resp.response_stages[-1].output_units = "count"
    resp.instrument_sensitivity.output_units = "count"

if not os.path.exists(OUTPUT_PATH):
    os.makedirs(OUTPUT_PATH)

# Begin response removal
files = sys.argv[1:]
print(f"removing response for {len(files)} waveforms")

# Check if filename in the right format for RENAME if chosen
if RENAME:
    try:
        net, sta, cha, loc, year, jul = os.path.basename(files[0]).split(".")
    except ValueError:
        print("error: RENAME=True but checked filename does not match format "
              f"NN.SSS.LL.CCC.YYYY.JJJ\n{os.path.basename(files[0])}")
        sys.exit(-1)

for fid in files:
    fid_out = os.path.basename(fid)
    print(fid_out, end="... ")

    # Check if this data has already been processed
    path_out = os.path.join(OUTPUT_PATH, fid_out)
    if os.path.exists(path_out) and not OVERWRITE:
        print("skipped, already processed")
        continue

    # Read data, ignore non waveform files
    try:
        st = read(fid)
    except (TypeError, IsADirectoryError):
        print("skipped, unknown file format")
        continue
  
    # Determine station naming from internal or from file
    if RENAME:
        net, sta, loc, cha, *_ = fid.split(".")
        # Rename internal stats based on filename
        st[0].stats.network = net
        st[0].stats.station = sta
        st[0].stats.location = loc
        st[0].stats.channel = cha
    else:
        net, sta, loc, cha = st[0].id.split(".")

    # Rename the inventory to match waveform so we can remove response
    inv[0].code = net
    inv[0][0].code = sta
    inv[0][0][0].code = cha
    inv[0][0][0].location_code = loc

    # Remove response with optional options
    st.remove_response(inventory=inv, pre_filt=PRE_FILT, output=OUTPUT)

    # Write out new file with response removed
    st.write(path_out, format="MSEED")
    print("done")

