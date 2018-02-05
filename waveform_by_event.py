"""5/2/18 Plot waveforms of temporary stations for a given event
"""
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from obspy import read
from obspy import read_inventory
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# event_id = False
# M4.8 near Wellington
# event_id = "2017p852531"
event_id = "2017p851921"
station_id = "RD06"
component = "HHZ"

# get event information
c = Client("GEONET")
if event_id:
    cat = c.get_events(eventid=event_id)
else:
    t_start = UTCDateTime("2017-320")
    t_end = UTCDateTime("2018-001")
    cat = c.get_events(starttime=t_start,
                        endtime=t_end,
                        minmagnitude=4,
                        maxmagnitude=6)

# get correct data filepaths
d1 = "/seis/prj/fwi/yoshi/RDF_Array/July2017_Sep2017/DATA_ALL/"
d2 = "/seis/prj/fwi/yoshi/RDF_Array/Sep2017_Nov2017/DATA_ALL/"
d3 = "/seis/prj/fwi/bchow/RDF_Array/Nov2017_Jan2018/DATA_ALL/"
# assuming event_id = True
origin = cat[0].origins[0].time
date_search = "{year}.{jday}".format(year=origin.year,jday=origin.julday)
filepath = []
for search in [d1,d2,d3]:
    filepath += glob.glob(search+"{sta}.{date}*{comp}*".format(sta=station_id,
                                                        date=date_search,
                                                        comp=component))

if not filepath:
    sys.exit("Files do doesn't not exist")

print("++{} data files found".format(len(filepath)))

resp_filepath = "/seis/prj/fwi/bchow/RDF_Array/DATALESS.RDF.XX"

# read in data - assuming only 1 datafile found
st = read(filepath[0])
inv = read_inventory(resp_filepath)

# preprocessing
# st.trim(starttime=origin-180,endtime=origin+3600*2)
st.taper(max_percentage=0.05)
st.attach_response(inv)
st.remove_response(output='DISP',water_level=200,plot=True)
st.detrend('simple')

# create time axis for plotting, initiate figure
t = np.linspace(0,st[0].stats.endtime-st[0].stats.starttime,st[0].stats.npts)
f,ax = plt.subplots(1)

# filter
for i in [3]:
    st_temp = st.copy()
    st_temp.filter('bandpass',freqmin=1/30,freqmax=1/i,corners=3,zerophase=True)
    plt.plot(t,st_temp[0].data,linewidth=1,label="{}s filter".format(i))

plt.legend()
plt.show()


import ipdb;ipdb.set_trace()
