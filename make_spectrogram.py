"""spectrogram of data for Kaikoura earthquake 2016-14-11T00:02:00
"""
import os
import sys
import glob
import obspy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from obspy.core.stream import Stream
from obspy import read
from obspy import read_inventory
from obspy import UTCDateTime

mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1

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
tmax = 50
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
f,(ax1,ax2,ax3) = plt.subplots(3,sharex=True,sharey=False,figsize=(9,5))
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
ax2.set_ylim([freqmin,0.2])
ax2.set_yticks(np.linspace(freqmin,0.2,3))
ax2.grid(color='w')

# subplot 2 (mirrored y axis for period) - this is impossible...
# locs,labels = plt.yticks()
# ax2a = ax2.twinx()
# ax2a.set_ylabel("Period (sec)")
# ax2a = ax2.yaxis.get_major_locator()
# # new_labels = 1/np.linspace(freqmin,0.2,10)
# # plt.yticks(locs,new_labels)
# # ax2a.set_yscale("log")
# import ipdb;ipdb.set_trace()

# processing for sub3 - determine duration of "secondary energy"
seismo =  st[0].data**2
samp_rate = int(st[0].stats.sampling_rate)
peak_amp = seismo.max()

# determine when 75% of energy reached
total_energy = np.trapz(seismo,t)
thresh_max = 0.9
thresh_min = 0.01
threshold = [thresh_min*total_energy, thresh_max*total_energy]

# iterate over waveform by 1 sec intervals
temp_energy = 0
energy_start = 0
while temp_energy < threshold[0]:
    energy_start += samp_rate
    temp_energy = np.trapz(seismo[0:energy_start],t[0:energy_start])

energy_end = energy_start
while temp_energy < threshold[1]:
    energy_end += samp_rate
    temp_energy = np.trapz(seismo[0:energy_end],t[0:energy_end])

duration = (energy_end-energy_start)/samp_rate

# subplot 3 edited seismogram with threshold
ax3.plot(t,seismo,'k',label=st[0].get_id())
ax3.plot(t[energy_start:energy_end],seismo[energy_start:energy_end],'r',
    label='{min}% < A^2 < {max}% total energy'.format(min=thresh_min*100,
                                                            max=thresh_max*100))

ax3.grid()
ax3.set_ylim([0,peak_amp+peak_amp*0.1])
ax3.set_ylabel('velocity^2 (m^2/s^2)')
ax3.set_xlabel("Time (sec)")
ax3.legend()

ano_x = round(t[int(energy_end + 100*samp_rate)]/1000) * 1000
ano_y = peak_amp/2
ax3.annotate("Duration of energy: {} sec".format(duration),
                                                xy=(ano_x,ano_y),
                                                xytext=(ano_x,ano_y))
# final touches
plt.xlim([200,2000])
plt.subplots_adjust(wspace=.5, hspace=0)

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
f.savefig(outpath,dpi=250)
print("saved figure:\n\t",outpath)

# text file containing parameters for easy comparisons
text_file = '/seis/prj/fwi/bchow/spectral/durations/{t0}-{t1}s_{p}p.txt'.format(
                                            t0=tmin,
                                            t1=tmax,
                                            p=int(thresh_max*100))
with open(text_file,'a') as f:
    write_string = ("\nID: {I}\n"
                    "Filter bounds: {F1}-{F2} s\n"
                    "Peak amplitude: {P} m/s\n"
                    "Duration: {D} s\n"
                    "Start: {S} s\n"
                    "End: {E} s\n".format(I=st[0].get_id(),
                                        F1=tmin,
                                        F2=tmax,
                                        P=st[0].data.max(),
                                        D=duration,
                                        S=energy_start/samp_rate,
                                        E=energy_end/samp_rate))
    f.write("="*80)
    f.write("{}".format(write_string))
    print("text file written")

# plt.show()
