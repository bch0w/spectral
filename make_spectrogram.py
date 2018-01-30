"""spectrogram of data for Kaikoura earthquake 2016-14-11T00:02:00
"""
import os
import sys
import glob
import obspy
import numpy as np
import matplotlib.pyplot as plt
from obspy.core.stream import Stream
from obspy import read
from obspy import read_inventory
from obspy import UTCDateTime

def geonet_data(station,comp,start,end=None,response=True):
    """
    returns a list of pathnames for GEONET archives on GNS internal system, if
    response == True, also returns path for response. If end not specified,
    returns a list of length 1 for the day in question

    :type station: str
    :param station: station name i.e. GKBS (case-insensitive)
    :type comp: str
    :param comp: component of interest, N/E/Z (case-insensitive)
    :type start: UTCDateTime
    :param start: starttime for data request
    :type end: UTCDateTime
    :param end: endtime for data request
    :type response: bool
    :param response: return response filepath
    """
    # channel naming convention based on instrument type
    chan_dict = {'Z':'HH','S':'BN'}
    channel = chan_dict[station[-1]] + chan_choice # i.e. HHZ

    # filepaths direct to GeoNet Archives path
    mseed_GNApath = '/geonet/seismic/{year}/NZ/{sta}/{ch}.D/'.format(year=year,
                                                                    sta=station,
                                                                    ch=channel)
    resp_GNApath = '/geonet/seed/RESPONSE/{}.NZ/'.format(station)

    # check if data spans multiple days
    if end:
        # if data spans multiple years

        # !!! DOESN'T WORK - need to figure out a way to dynamically determine
        #!!! which filenames to grab. Regex maybe? want this to be a general fx
        if start.year != end.year:
            mseed_files = []
            for iterate_year in range(start.year,end_year+1,1):
                mseed_files += glob.glob(mseed_GNApath + '*{year}.*'.format(
                                                            year=start.year))
            mseed_files.sort()
    else:
        mseed_files = glob.glob(mseed_GNApath +'*{year}.{day}'.format(
                                                            year=start.year,
                                                            day=start.julday))

    # mseed sets naming parameters for response
    NET,STA,LOC,CHA,SUFX,YEAR,JDAY = os.path.basename(
                                                mseed_files[0]).split('.')

    if response:
        response_filename = "RESP.{net}.{sta}.{loc}.{cha}".format(net=NET,
                                                                    sta=STA,
                                                                    loc=LOC,
                                                                    cha=CHA)
        response_filepath = os.path.join(resp_GNApath,response_filename)


# ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# set station and channel ID
station = sys.argv[1].upper() # i.e. gkbs
chan_choice = sys.argv[2].upper() # i.e. z
if chan_choice not in ['N','E','Z','1','2']:
    sys.exit('incorrect channel choice')

chan_dict = {'Z':'HH','S':'BN'}
channel = chan_dict[station[-1]] + chan_choice # i.e. HHZ

# origin time
kaikoura_origin_time = UTCDateTime('2016-11-13T11:02:56')
start_time = kaikoura_origin_time - 200
end_time = kaikoura_origin_time + 2000

year = start_time.year
julday = start_time.julday

# set filepaths direct to geonet archives path (GNA)
# base_filepath = '/seis/prj/fwi/bchow/geonet_stations/NZ/' # old
mseed_GNApath = '/geonet/seismic/{year}/NZ/{sta}/{ch}.D/'.format(year=year,
                                                            sta=station,
                                                            ch=channel)
resp_GNApath = '/geonet/seed/RESPONSE/{}.NZ/'.format(station)

# check if data spans multiple days
if start_time.julday != julday:
    multiday = True
    mseed_filename = glob.glob(mseed_GNApath +'*{}'.format(julday-1))
    mseed_filename += glob.glob(mseed_GNApath +'*{}'.format(julday))
    mseed_filename.sort()
else:
    multiday = False
    mseed_filename = glob.glob(mseed_GNApath +'*{}'.format(julday))

# read in response (filename set naming parameters)
NET,STA,LOC,CHA,SUFX,YEAR,JDAY = os.path.basename(mseed_filename[0]).split('.')
response_filename = "RESP.{net}.{sta}.{loc}.{cha}".format(net=NET,
                                                            sta=STA,
                                                            loc=LOC,
                                                            cha=CHA)
response_filepath = os.path.join(resp_GNApath,response_filename)
inv = read_inventory(response_filepath)

# read in data streams
if multiday:
    tr1 = read(mseed_filename[0])[0]
    tr2 = read(mseed_filename[1])[0]
    st = Stream(tr1+tr2)
    st_original = st.copy()
else:
    st = read(mseed_filename[0])
    st_original = st.copy()

# trim
st.trim(starttime=start_time,endtime=end_time)

# remove response
st.attach_response(inventories=inv)
st.remove_response(output='VEL',water_level=100)

# filter window
tmin = 1
tmax = 100
freqmin = 1/tmax
freqmax = 1/tmin

# spectrogram settings
window_length = 200
clip_low = 0
clip_high = 1.0

# preprocessing
st.detrend("simple")
st.filter("bandpass",freqmin=freqmin,freqmax=freqmax,corners=3)
st.taper(max_percentage=0.05)

# plotting
print("plotting",end=", ")
f,(ax1,ax2) = plt.subplots(2,sharex=True,sharey=False,figsize=(9,5))
label_dict = {'Z':'Vertical','1':'N/S','N':'N/S','2':'E/W','E':'E/W'}

# subplot 1 (waveform)
t = np.linspace(0,st[0].stats.endtime-st[0].stats.starttime,st[0].stats.npts)
ax1.plot(t,st[0].data,linewidth=1,c='k',label=st[0].get_id())
ax1.legend()
ax1.grid(which="both")
ax1.set_ylabel("{} Velocity (m/s)".format(label_dict[st[0].stats.channel[-1]]))
plot_title = "Kaikoura Earthquake | " + \
"{instr} | {otime} | Filt: {tmin}-{tmax} s | FFT Len: {wl}s | Clip: {c0},{c1}".format(
            otime=start_time,
            instr=st[0].get_id(),
            tmin=tmin,
            tmax=tmax,
            wl=window_length,
            c0=clip_low,
            c1=clip_high
            )
ax1.set_title(plot_title,fontsize=10)

# subplot 2 (spectrogram)
spec = st.spectrogram(log=False,clip=[clip_low,clip_high],cmap='plasma',
                        axes=ax2,show=False,wlen=window_length)
ax2.set_ylabel("Frequency (Hz)")
ax2.set_xlabel("Time (sec)")
ax2.set_ylim([freqmin,0.2])
ax2.set_yticks(np.linspace(freqmin,0.2,10),10)
ax2.grid(color='w')

# subplot 2 (mirrored y axis for period) - Didn't work properly
# ax2a = ax2.twinx()
# ax2a.set_ylabel("Period (sec)")
# ax2a.set_yticks(np.linspace(tmin,tmax,10),10)
# ax2a.set_ylim([tmax,tmin])
# ax2a.set_yscale("log")

# final touches
plt.xlim([200,2000])
plt.subplots_adjust(wspace=0, hspace=0)

# save fig
figure_folder = '/seis/prj/fwi/bchow/spectral/figures/spectrograms/'
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
# f.savefig(outpath,dpi=f.dpi)
plt.show()
print("saved figure:\n\t",outpath)

# plt.show()
