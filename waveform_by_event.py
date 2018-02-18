"""5/2/18 Plot waveforms of temporary stations for a given event
"""
import os
import sys
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
from getdata import geonet_internal, pathnames, vog, fdsn_download
from obspy import read, read_inventory, UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import locations2degrees

# ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# user input arguments
parser = argparse.ArgumentParser(description='Plot waveforms for earthquake\
data, example call $ python waveform_by_event.py --station TBAS --channel\
 BNZ --end 2015-02-01 --unit vel')
parser.add_argument('--choice', help='GEONET or RDF, default = RDF',
                    type=str,default='RDF')
parser.add_argument('--station', help='Station choice, default = RD01',type=str,
                    default='RD01')
parser.add_argument('--channel', help='Instrument channel for choice geonet\
                    , i.e. BN1/HHZ, ignores component, default = HHZ',
                    type=str, default='HHZ')
parser.add_argument('--unit', help='Units of waveforms, [vel]/disp/acc',
                    type=str,default='VEL')
parser.add_argument('--start', help='Starttime for catalog, default=\
                    2017-320',type=str, default='2017-320')
parser.add_argument('--end', help='Endtime, default = 2018-001',type=str,
                    default='2018-001')

# parse arguments
arg = parser.parse_args()
choice = arg.choice.lower() # geonet or rdf
station_id = arg.station.upper() # i.e. RD06 or PXZ
channel = arg.channel.upper()
response_output = arg.unit.upper() # i.e. "VEL" or "DISP" or "ACC"
t_start = UTCDateTime(arg.start)
t_end = UTCDateTime(arg.end)

vic_or_gns = vog()

# list of earthquakes to look at
event_list = ["2017p852531", # M4.8 36km; W of Wellington
              "2017p851921", # M4.1 39km; NW of Palmy, might be too close
              "2017p860319", # M4.3 15km; NE of Masterton
              "2017p861155", # M4.7 11km; SW of Wellington
              "2017p968142", # M4.1 16km; W of Mt Ruapehu
              "2017p879247"] # M4.7 40km; N of Bay of Plenty
# event_id = "2017p916322"
event_id = False

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

print("{} events found".format(len(cat)))
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

        st = read(path_vert[0]) + read(path_north[0]) + read(path_east[0])
        resp_filepath = pathnames(vic_or_gns)['rdf'] + "DATALESS.RDF.XX"
        inv = read_inventory(resp_filepath)

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

    # grab geonet data either internally or via fdsn
    elif choice == "geonet":
        if vic_or_gns == "GNS":
            path_vert, resp_vert = geonet_internal(station=station_id,
                                        channel= channel[:2] + 'Z',
                                        start = origin,
                                        response = True)
            path_north, resp_north = geonet_internal(station=station_id,
                                        channel= channel[:2] + 'N',
                                        start = origin,
                                        response = True)
            path_east, resp_east = geonet_internal(station=station_id,
                                        channel= channel[:2] + 'E',
                                        start = origin,
                                        response = True)
            inv = read_inventory(resp_vert)
            inv += read_inventory(resp_north)
            inv += read_inventory(resp_east)
            st = read(path_vert[0]) + read(path_north[0]) + read(path_east[0])


        elif vic_or_gns == "VIC":
            st, inv = fdsn_download(station=station_id,
                                                channel = channel[:2] + '*',
                                                start = origin,
                                                response = True)

    # preprocessing
    pushback = -50
    # st.merge()
    st.trim(starttime=origin+pushback,endtime=origin+60*7.5)
    st.detrend('linear')
    st.taper(max_percentage=0.05)
    st.attach_response(inv)
    st.remove_response(output=response_output,water_level=60)

    # create time axis for plotting, initiate figure
    try:
        stats = st[0].stats
    except IndexError:
        # if data massing, skip over
        print(origin)
        continue
    t = np.linspace(pushback,stats.endtime-stats.starttime,stats.npts)
    f,(ax1,ax2,ax3) = plt.subplots(3,figsize=(9,5),
                                     sharex=True,
                                     sharey=True,
                                     dpi=200)

    # clean up trace
    st_temp = st.copy()
    st_temp.detrend('linear')
    st_temp.taper(max_percentage=0.05)
    if channel[0] == 'E':
        hp_cutoff = 1/10
    else:
        hp_cutoff = 1/30
    st_temp.filter('highpass',freq=hp_cutoff,corners=3,zerophase=True)

    # signal to noise ratio
    samp_rate = int(st_temp[0].stats.sampling_rate)
    noise_data = st_temp.select(component='Z')[0].data[0:-pushback*samp_rate]
    noise_data = np.sqrt(noise_data**2)
    noise_mean = noise_data.mean()
    peak_signal = st_temp.select(component='Z')[0].data.max()
    signal_noise_ratio = peak_signal/noise_mean
    print("Signal to noise ratio ",signal_noise_ratio)

    # filter bands
    alpha = 0.4
    linewidth = 0.5
    for i in [1,3,6,10]:
        for ax,co in zip([ax1,ax2,ax3],['N','E','Z']):

            # for each filter band, plot each component trace
            st_filter = st_temp.copy()
            st_filter.filter('lowpass',freq=1/i,corners=3,zerophase=True)

            # plot the raw seismograms in the background, only once
            if i == 1:
                ax.plot(t,st_temp.select(component=co)[0].data,
                            linewidth=0.4,
                            alpha=0.2,
                            label="raw".format(round(1/i,2),i),
                            zorder=i,
                            color='gray')

            # plot filtered waveform
            ax.plot(t,st_filter.select(component=co)[0].data,
                        linewidth=linewidth,
                        alpha=alpha,
                        label="{}Hz/{}s lowpass".format(round(1/i,2),i),
                        zorder=i+1)
            if i == 3:
                ax.set_ylim([min(st_filter[0].data),max(st_filter[0].data)])


        alpha+=0.2
        linewidth+=0.2
        # print("Lowpass at {} s".format(st_filter[0].data.max()))

    # ax.set_ylim([-0.000023,0.000023])
    # determine ylimits with min/max raw trace values
        # maxs,mins = [],[]
        # for co in ['N','E','Z']:
        #     maxs.append(max(st_temp.select(component=co)[0].data))
        #     mins.append(min(st_temp.select(component=co)[0].data))
        # ymin = min(mins)
        # ymax= max(maxs)
        # ax.set_ylim([ymin,ymax])
        # print("min: {} counts/max: {} counts".format(ymin,ymax))


    ax3.set_xlabel('Time (sec)')
    unit_dict = {"VEL":"m/s","DISP":"m"}
    ax1.set_ylabel('North {} ({})'.format(response_output,
                                            unit_dict[response_output]))
    ax2.set_ylabel('East {} ({})'.format(response_output,
                                            unit_dict[response_output]))
    ax3.set_ylabel('Vertical {} ({})'.format(response_output,
                                            unit_dict[response_output]))
    ax1.set_title("{s} | {e} | {o} | M{m}".format(s=station_id,
                                                            e=event_id,
                                                            o=origin,
                                                            m=magnitude))
                                                            # d=distance))
    for ax in [ax1,ax2,ax3]: ax.grid(True)
    ax2.legend(loc='best',prop={'size':7})
    ax2.set_zorder(100)
    basepath = pathnames(vic_or_gns)['plots'] + "waveforms/"
    if not os.path.exists(basepath+event_id):
        os.makedirs(basepath+event_id)
    figure_name = os.path.join(basepath,event_id,"{0}_{1}_{2}.png".format(
                                                event_id,
                                                station_id,
                                                response_output))
    # print(figure_name)
    # plt.subplots_adjust(wspace=.5, hspace=0)
    # plt.savefig(figure_name)

    plt.show()
