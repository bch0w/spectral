"""5/2/18 Plot waveforms of temporary stations for a given event
"""
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from getdata import geonet_data
from obspy import read
from obspy import read_inventory
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# list of earthquakes to look at
event_list = ["2017p852531", # M4.8 36km; W of Wellington
              "2017p851921", # M4.1 39km; NW of Palmy, might be too close
              "2017p860319", # M4.3 15km; NE of Masterton
              "2017p861155", # M4.7 11km; SW of Wellington
              "2017p968142", # M4.1 16km; W of Mt Ruapehu
              "2017p879247"] # M4.7 40km; N of Bay of Plenty
event_id = event_list[4]

print(event_id)

try:
    choice = sys.argv[1].lower() # geonet or temp
    station_id = sys.argv[2].upper() # i.e. RD06 or PXZ
    component = sys.argv[3].upper() # i.e. HHZ
    response_output = sys.argv[4].upper() # i.e. "VEL" or "DISP"
except IndexError:
    sys.exit("example call: python waveform_by_event.py temp RD06 HHZ VEL")

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

# assuming event_id is True
origin = cat[0].origins[0].time
# origin = UTCDateTime('2017-11-12:18:03:58')

# temporary station filepaths
if choice == "temp":
    d1 = "/seis/prj/fwi/yoshi/RDF_Array/July2017_Sep2017/DATA_ALL/"
    d2 = "/seis/prj/fwi/yoshi/RDF_Array/Sep2017_Nov2017/DATA_ALL/"
    d3 = "/seis/prj/fwi/bchow/RDF_Array/Nov2017_Jan2018/DATA_ALL/"
    date_search = "{year}.{jday}".format(year=origin.year,jday=origin.julday)
    filepath = []
    for search in [d1,d2,d3]:
        filepath += glob.glob(search+"{sta}.{date}*{comp}*".format(
                                                            sta=station_id,
                                                            date=date_search,
                                                            comp=component))

    if not filepath:
        sys.exit("Files do doesn't not exist")

    print("++ {} data files found".format(len(filepath)))
    resp_filepath = "/seis/prj/fwi/bchow/RDF_Array/DATALESS.RDF.XX"

elif choice == "geonet":
    filepath,resp_filepath = geonet_data(station=station_id,
                                        comp=component[-1],
                                        start=origin)


# read in data - assuming only 1 datafile found
st = read(filepath[0])
inv = read_inventory(resp_filepath)

# preprocessing
st.trim(starttime=origin-180,endtime=origin+60*5)
st.detrend('linear')
st.taper(max_percentage=0.05)
st.attach_response(inv)
pre_filt = [0.025,0.033,1,10]
st.remove_response(output=response_output,water_level=60)

# create time axis for plotting, initiate figure
t = np.linspace(0,st[0].stats.endtime-st[0].stats.starttime,st[0].stats.npts)
f,ax = plt.subplots(1,figsize=(9,5))

# clean up trace
st_temp = st.copy()
st_temp.detrend('linear')
st_temp.taper(max_percentage=0.05)
st_temp.filter('highpass',freq=1/30)

# no lowpass
st_nolowpass = st_temp.copy()
plt.plot(t,st_nolowpass[0].data,linewidth=0.5,label="No lowpass")
# print("No lowpass: {}".format(st_nolowpass[0].data.max()))

# filter bands
for i in [3,6,10]:
    st_filter = st_temp.copy()
    st_filter.filter('lowpass',freq=1/i)
    plt.plot(t,st_filter[0].data,
                linewidth=0.5,
                label="{}s filter".format(i),
                zorder=i)
    print("Lowpass at {} s".format(st_filter[0].data.max()))

plt.xlabel('Time (sec)')
unit_dict = {"VEL":"m/s","DISP":"m"}
plt.ylabel('{} ({})'.format(response_output,unit_dict[response_output]))
plt.title("Station {s} | Event {e} | Origin {o}".format(s=station_id,
                                                        e=event_id,
                                                        o=origin))
plt.grid(True)
plt.legend()
basepath = "/seis/prj/fwi/bchow/spectral/output_plots/waveforms/"
if not os.path.exists(basepath+event_id):
    os.makedirs(basepath+event_id)
figure_name = os.path.join(basepath,event_id,"{0}_{1}.{2}_{3}.png".format(
                                            event_id,
                                            station_id,
                                            component,
                                            response_output))
print(figure_name)
plt.savefig(figure_name)
# plt.show()
