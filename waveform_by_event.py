"""5/2/18 Plot waveforms of temporary stations for a given event
"""
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from getdata import geonet_data, pathnames
from obspy import read, read_inventory, UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import locations2degrees

vic_or_gns = 'vic'

# ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# does this even work?
# plt.rcParams['axes.axisbelow'] = True

# list of earthquakes to look at
event_list = ["2017p852531", # M4.8 36km; W of Wellington
              "2017p851921", # M4.1 39km; NW of Palmy, might be too close
              "2017p860319", # M4.3 15km; NE of Masterton
              "2017p861155", # M4.7 11km; SW of Wellington
              "2017p968142", # M4.1 16km; W of Mt Ruapehu
              "2017p879247"] # M4.7 40km; N of Bay of Plenty
event_id = False

try:
    choice = sys.argv[1].lower() # geonet or rdf
    station_id = sys.argv[2].upper() # i.e. RD06 or PXZ
    response_output = sys.argv[3].upper() # i.e. "VEL" or "DISP"
except IndexError:
    sys.exit("python waveform_by_event.py rdf/geonet RD??/PXZ VEL/DISP")

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
                        maxmagnitude=6,
                        minlatitude=-50,
                        maxlatitude=-35,
                        minlongitude=165,
                        maxlongitude=180,
                        orderby="magnitude")

for event in cat:
    # grab earthquake data - magnitude is M, ignore Ml and Mlv
    event_id = event.resource_id.id.split('/')[1]
    origin = event.origins[0].time
    magnitude = event.magnitudes[2].mag
    latitude = event.origins[0].latitude
    longitude = event.origins[0].longitude

    # temporary station filepaths
    if choice == "rdf":
        d1 = pathnames(vic_or_gns)['rdf'] + "July2017_Sep2017/DATA_ALL/"
        d2 = pathnames(vic_or_gns)['rdf'] + "Sep2017_Nov2017/DATA_ALL/"
        d3 = pathnames(vic_or_gns)['rdf'] + "Nov2017_Jan2018/DATA_ALL/"
        date_search = "{year}.{jday}".format(
                                            year=origin.year,jday=origin.julday)
        path_vert,path_north,path_east = [],[],[]
        for search in [d1,d2,d3]:
            path_vert += glob.glob(search+"{sta}.{date}*HHZ*".format(
                                                            sta=station_id,
                                                            date=date_search))
            path_north += glob.glob(search+"{sta}.{date}*HHN*".format(
                                                            sta=station_id,
                                                            date=date_search))
            path_east += glob.glob(search+"{sta}.{date}*HHE*".format(
                                                            sta=station_id,
                                                            date=date_search))

        # print("++ {} data files found".format(len(filepath)))
        resp_filepath = pathnames(vic_or_gns)['rdf'] + "DATALESS.RDF.XX"

    elif choice == "geonet":
        filepath,resp_filepath = geonet_data(station=station_id,
                                            comp=component[-1],
                                            start=origin)

    # read in data - assuming only 1 datafile found
    st = read(path_vert[0]) + read(path_north[0]) + read(path_east[0])
    # incase data gaps (??? is this okay?)
    st.merge()
    # inv = read_inventory(resp_filepath)

    # grab station information - determine station-event distance
    stations,lats,lons = [],[],[]
    with open(pathnames(vic_or_gns)['rdf']+ 'rdf_locations.txt','r') as f:
        station_lines = f.readlines()
    for lines in station_lines[1:]:
        splitlines = lines.split(',')
        if splitlines[0] == station_id:
            distance = locations2degrees(lat1=latitude,
                                        long1=longitude,
                                        lat2=float(splitlines[2]),
                                        long2=float(splitlines[3]))
            distance *= 40000/360
            distance = round(distance,2)

    # preprocessing
    pushback = -30
    st.trim(starttime=origin+pushback,endtime=origin+60*7.5)
    st.detrend('linear')
    st.taper(max_percentage=0.05)
    # st.attach_response(inv)
    # st.remove_response(output=response_output,water_level=60)

    # create time axis for plotting, initiate figure
    stats = st[0].stats
    t = np.linspace(pushback,stats.endtime-stats.starttime,stats.npts)
    f,(ax1,ax2,ax3) = plt.subplots(3,figsize=(9,5),
                                     sharex=True,
                                     sharey=True,
                                     dpi=200)

    # clean up trace
    st_temp = st.copy()
    st_temp.detrend('linear')
    st_temp.taper(max_percentage=0.05)
    st_temp.filter('highpass',freq=1/30,corners=3,zerophase=True)

    # filter bands
    alpha = 0.7
    linewidth = 0.5
    for i in [1,3,6,10]:
        for ax,co in zip([ax1,ax2,ax3],['N','E','Z']):
            # for each filter band, plot each component trace
            # st_filter = st_temp.copy()
            # st_filter.filter('lowpass',freq=1/i,corners=3,zerophase=True)
            # ax.plot(t,st_filter.select(component=co)[0].data,
            #             linewidth=linewidth,
            #             alpha=alpha,
            #             label="{}Hz/{}s lowpass".format(round(1/i,2),i),
            #             zorder=i)
            # plot the raw seismograms in the background, only once
            if i == 1:
                ax.plot(t,st_temp.select(component=co)[0].data,
                            linewidth=0.4,
                            # alpha=0.2,
                            label="raw".format(round(1/i,2),i),
                            zorder=100,
                            color='k')
                # ax.set_ylim([max(st_filter[tr].data),min(st_filter[tr].data)])

        alpha+=0.1
        linewidth+=0.2
        # print("Lowpass at {} s".format(st_filter[0].data.max()))

    # determine ylimits with min/max trace values
    maxs,mins = [],[]
    for co in ['N','E','Z']:
        maxs.append(max(st_temp.select(component=co)[0].data))
        mins.append(min(st_temp.select(component=co)[0].data))
    ymin = min(mins)
    ymax= max(maxs)
    ax.set_ylim([ymin,ymax])
    print("min: {} counts/max: {} counts".format(ymin,ymax))


    ax3.set_xlabel('Time (sec)')
    unit_dict = {"VEL":"m/s","DISP":"m"}
    ax1.set_ylabel('North {} ({})'.format(response_output,
                                            unit_dict[response_output]))
    ax2.set_ylabel('East {} ({})'.format(response_output,
                                            unit_dict[response_output]))
    ax3.set_ylabel('Vertical {} ({})'.format(response_output,
                                            unit_dict[response_output]))
    ax1.set_title("{s} | {e} | {o} | M{m} | ~{d}km".format(s=station_id,
                                                            e=event_id,
                                                            o=origin,
                                                            m=magnitude,
                                                            d=distance))
    for ax in [ax1,ax2,ax3]: ax.grid(True)
    ax2.legend(loc='best',prop={'size':7})
    basepath = pathnames(vic_or_gns)['plots'] + "waveforms/"
    if not os.path.exists(basepath+event_id):
        os.makedirs(basepath+event_id)
    figure_name = os.path.join(basepath,event_id,"{0}_{1}_{2}.png".format(
                                                event_id,
                                                station_id,
                                                response_output))
    # print(figure_name)
    # plt.savefig(figure_name)
    plt.subplots_adjust(wspace=.5, hspace=0)

    plt.show()
