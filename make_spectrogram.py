"""spectrogram of data for Kaikoura earthquake 2016-14-11T00:02:00
only works on GNS computer
"""
import os
import sys
import glob
import obspy
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from obspy.core.stream import Stream
from obspy import read, read_inventory, UTCDateTime
from obspy.clients.fdsn import Client
from getdata import vog, geonet_internal, fdsn_download, pathnames


mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1
# ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1

def event_stream(vic_or_gns,choice,station,channel,event_id=False):
    """Given a GEONET event ID, return waveform streams with removed
    response. Waveforms can be from RDF or GEONET permanent stations, chooses
    correct downloading format based on requests
    :type vic_or_gns: str
    :param vic_or_gns: specifies pathames for data searching
    :type choice: str
    :param choice: GEONET or RDF (temporary) station data
    :type station: str
    :param station: station name i.e. GKBS or RD01 (case-insensitive)
    :type channel: str
    :param channel: channel of interest, i.e. BHN, HHE (case-insensitive)
    :type event_id: str
    :param event_id: GEONET event id for earthquakes, if default, function will
                    simply collect a catalog of events for a given time periods
    """
    # parse arguments
    choice = choice.lower() # geonet or rdf
    station_id = station.upper() # i.e. RD06 or PXZ
    channel = channel.upper()

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
        origin = event.origins[0].time

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
                st = (read(path_vert[0]) + read(path_north[0]) +
                                                            read(path_east[0]))


            elif vic_or_gns == "VIC":
                st, inv = fdsn_download(station=station_id,
                                        channel = channel[:2] + '*',
                                        start = origin,
                                        response = True)

        return st, inv, cat


# MAIN
parser = argparse.ArgumentParser(description='Plot spectrograms for earthquake\
                    example call python make_spectrogram.py --choice rdf\
                     --station rd01 --channel HHZ --output vel')
parser.add_argument('--choice', help='GEONET or RDF, default = RDF',
                    type=str,default='RDF')
parser.add_argument('--station', help='Station choice, default = RD01',type=str,
                    default='RD01')
parser.add_argument('--channel', help='Instrument channel for choice geonet\
                    , i.e. BN1/HHZ, ignores component, default = HHZ',
                    type=str, default='HHZ')
parser.add_argument('--output', help='Units of waveforms, [vel]/disp/acc',
                    type=str,default='VEL')
parser.add_argument('--component', help='Component of waveforms, N/E/Z',
                    type=str,default='Z')
parser.add_argument('--event_id', help='Quake event id from GEONET',
                    type=str, default='2017p916322')


# parse arguments
arg = parser.parse_args()
choice = arg.choice.lower() # geonet or rdf
station_id = arg.station.upper() # i.e. RD06 or PXZ
channel = arg.channel.upper()
output = arg.output.upper() # i.e. "VEL" or "DISP" or "ACC"
comp = arg.component.upper()

# grab data
vic_or_gns = vog()
st,inv,cat = event_stream(vic_or_gns=vic_or_gns,
                      choice=choice,
                      station=station_id,
                      channel=channel,
                      event_id="2017p916322")

# origin time
event = cat[0]
event_id = event.resource_id.id.split('/')[1]
origin = event.origins[0].time
start_time = origin - 200
end_time = origin + 600

# filter window
tmin = 3
tmax = 30
freqmin = 1/tmax
freqmax = 1/tmin

# spectrogram settings
window_length = 200
clip_low = 0
clip_high = 1.0

# preprocessing
st.detrend("simple")
st.taper(max_percentage=0.05)
st.filter("bandpass",freqmin=freqmin,freqmax=freqmax,corners=3)

# separate horizontal data streams
north = st.select(component='1')
east = st.select(component='2')
if (len(north) and len(east)) == 0:
    north = st.select(component='N')
    east = st.select(component='E')
north = north[0].data
east = east[0].data

# choices for working trace
vertical = st.select(component='Z')[0].data
horizontal = np.sqrt(north**2 + east**2)
motion_vector = np.sqrt(vertical**2 + north**2 + east**2)

# assign trace to use
tr_work = st[0].copy()
tr_work.data = motion_vector
stats = tr_work.stats

# plotting
f,(ax1,ax1a,ax1b,ax2,ax3) = plt.subplots(5,
                                        sharex=True,
                                        sharey=False,
                                        figsize=(9,5))

# s============= subplot 1 (waveform) =============
tr = st.select(component=comp)[0]
stats  = tr.stats

t = np.linspace(0,stats.endtime-stats.starttime,stats.npts)
ax1.plot(t,tr.data,linewidth=1,c='k',label=tr.get_id())
ax1.legend()
ax1.grid(which="both")
ax1.set_ylabel("{} Velocity (m/s)".format(label_dict[stats.channel[-1]]))
plot_title = "{id} | {i} | {o} | {t0}-{t1} s | FFT Len: {wl}s".format(
            id=event_id,
            i=tr.get_id(),
            o=start_time,
            t0=tmin,
            t1=tmax,
            wl=window_length
            )
ax1.set_title(plot_title,fontsize=10)

# ============= subplot 2 (spectrogram) =============
spec = tr.spectrogram(log=False,clip=[clip_low,clip_high],cmap='plasma',
                        axes=ax2,show=False,wlen=window_length)
ax2.set_ylabel("Frequency (Hz)")
ax2.set_ylim([freqmin,2])
# ax2.set_yticks(np.linspace(freqmin,0.2,3))
ax2.grid(color='w')


# ===== processing for sub3 - determine duration of "secondary energy" ======
seismo =  tr.data**2
samp_rate = stats.sampling_rate

peak_amp = seismo.max()
threshold_percentage = 0.1
threshold = peak_amp * threshold_percentage

# loop over seismogram, determine start and end of peak energy
a_over, s_over = [],[]
for i,amplitude in enumerate(seismo):
    if amplitude >= threshold:
        s_over.append(i)
        a_over.append(amplitude)

# find edgepoints by checking if the next sample j is the same as i+1
s_edge,a_edge,sections = [s_over[0]],[a_over[0]],[]
for i,(S,A) in enumerate(zip(s_over[1:-2],a_over[1:-2])):
    if s_over[i+2] != (S + 1):
        section = np.trapz(a_edge,s_edge)
        sections.append(section)
        s_edge,a_edge = [s_over[i+1]],[a_over[i+1]]
    else:
        s_edge.append(S)
        a_edge.append(A)

# convert samples to time
t_over = []
for S in s_over:
    t_over.append(t[S])

duration = sum(sections)
print("Duration criteria: {}".format(duration))

# subplot 3 edited seismogram with threshold
ax3.plot(t,seismo,'k',label=st[0].get_id())
ax3.scatter(t_over,a_over,c='r',marker='.',s=1,zorder=100)
h_lab = "Threshold = {}% peak amplitude".format(int(threshold_percentage*100))
ax3.axhline(y=threshold,
            xmin=t[0],xmax=t[-1],
            zorder=1,
            color='r',
            linestyle='-.',
            linewidth=.85,
            label=h_lab)

ax3.grid()
ax3.set_ylim([0,peak_amp+peak_amp*0.1])
ax3.set_ylabel('velocity^2 (m^/s^2)')
ax3.set_xlabel("Time (sec)")
ax3.legend()

# ano_x = round(t[int(energy_end + 100*samp_rate)]/1000) * 1000
# ano_y = threshold * 2
# ax3.annotate("Duration of energy: {} sec".format(duration),
#                                                 xy=(ano_x,ano_y),
#                                                 xytext=(ano_x,ano_y))

# final touches
plt.xlim([200,2000])
plt.subplots_adjust(wspace=.5, hspace=0)

# save fig
figure_folder = '/seis/prj/fwi/bchow/spectral/output_plots/spectrograms/'
subfolder = '{}-{}s'.format(tmin,tmax)
foldercheck = os.path.join(figure_folder,subfolder)
if not os.path.exists(foldercheck):
    print("Making directory ",foldercheck,end=", ")
    os.makedirs(foldercheck)
figure_name = "{s}-{c}_{l}-{h}kaikoura".format(s=station,
                                c=channel,
                                l=tmin,
                                h=tmax)
outpath = os.path.join(figure_folder,subfolder,figure_name)
# f.savefig(outpath,dpi=250)
# print("saved figure:\n\t",outpath)

# text file containing parameters for easy comparisons
text_file = '/seis/prj/fwi/bchow/spectral/duration/{t0}-{t1}s_amp.txt'.format(
                                            t0=tmin,
                                            t1=tmax)
# with open(text_file,'a') as f:
#     write_string = ("\nID: {I}\n"
#                     "Filter bounds: {F1}-{F2} s\n"
#                     "Peak amplitude: {P} m/s\n"
#                     "Duration: {D} s\n"
#                     "Start: {S} s\n"
#                     "End: {E} s\n".format(I=st[0].get_id(),
#                                         F1=tmin,
#                                         F2=tmax,
#                                         P=st[0].data.max(),
#                                         D=duration,
#                                         S=energy_start/samp_rate,
#                                         E=energy_end/samp_rate))
#     f.write("="*80)
#     f.write("{}".format(write_string))
#     print("text file written")

plt.show()
