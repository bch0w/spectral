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
st.trim(starttime=start_time,endtime=end_time)
st.detrend('linear')
st.taper(max_percentage=0.05)
st.attach_response(inventories=inv)
st.remove_response(output=output,water_level=100)

st.detrend('linear')
st.taper(max_percentage=0.05)
st.filter("bandpass",freqmin=freqmin,freqmax=freqmax,corners=3)

# plotting
print("plotting",end=", ")
f,(ax1,ax2,ax3) = plt.subplots(3,sharex=True,sharey=False,figsize=(9,5))
label_dict = {'Z':'Vertical','1':'N/S','N':'N/S','2':'E/W','E':'E/W'}

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
threshold_percentage = 0.05
threshold = peak_amp * threshold_percentage

# loop over seismogram, determine start and end of peak energy
over_threshold = []
for i,sample in enumerate(seismo):
    if sample >= threshold:
        over_threshold.append(i)

energy_start = over_threshold[0]
energy_end = over_threshold[-1]
duration = (energy_end-energy_start)/samp_rate

# ============= subplot 3 edited seismogram with threshold =============
ax3.plot(t,seismo,'k',label=tr.get_id())
ax3.plot(t[energy_start:energy_end],seismo[energy_start:energy_end],'r',
                                                label='Amplitude >= Threshold')
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

ano_x = round(t[int(energy_end + 100*samp_rate)]/1000) * 1000
ano_y = threshold * 2
ax3.annotate("Duration of energy: {} sec".format(duration),
                                                xy=(ano_x,ano_y),
                                                xytext=(ano_x,ano_y))


# ============= final touches =============
# plt.xlim([200,2000])
plt.subplots_adjust(wspace=.5, hspace=0)

plt.show()


# ===========================OLD-INTEGRAL-APPROACH==============================
    # # determine when 75% of energy reached
    # total_energy = np.trapz(seismo,t)
    # thresh_max = 0.9
    # thresh_min = 0.01
    # threshold = [thresh_min*total_energy, thresh_max*total_energy]
    #
    # # iterate over waveform by 1 sec intervals
    # temp_energy = 0
    # energy_start = 0
    # while temp_energy < threshold[0]:
    #     energy_start += samp_rate
    #     temp_energy = np.trapz(seismo[0:energy_start],t[0:energy_start])
    #
    # energy_end = energy_start
    # while temp_energy < threshold[1]:
    #     energy_end += samp_rate
    #     temp_energy = np.trapz(seismo[0:energy_end],t[0:energy_end])
    #
    # duration = (energy_end-energy_start)/samp_rate
    #
    # # subplot 3 edited seismogram with threshold
    # ax3.plot(t,seismo,'k',label=st[0].get_id())
    # ax3.plot(t[energy_start:energy_end],seismo[energy_start:energy_end],'r',
    #     label='{min}% < A^2 < {max}% total energy'.format(min=thresh_min*100,
    #                                                             max=thresh_max*100))
    #
    # ax3.grid()
    # ax3.set_ylim([0,peak_amp+peak_amp*0.1])
    # ax3.set_ylabel('velocity^2 (m^2/s^2)')
    # ax3.set_xlabel("Time (sec)")
    # ax3.legend()
    #
    # ano_x = round(t[int(energy_end + 100*samp_rate)]/1000) * 1000
    # ano_y = peak_amp/2
    # ax3.annotate("Duration of energy: {} sec".format(duration),
    #                                                 xy=(ano_x,ano_y),
    #                                                 xytext=(ano_x,ano_y))
# ===========================OLD-INTEGRAL-APPROACH==============================
