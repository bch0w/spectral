"""spectrogram of data for Kaikoura earthquake 2016-14-11T00:02:00
only works on GNS computer
"""
import os
import sys
import glob
import obspy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from obspy.core.stream import Stream
from obspy import read, read_inventory, UTCDateTime
from getdata import geonet_internal

# ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1

# set station and channel ID
channel_dict = {"Z":"HH*","S":"BN*"}
station = sys.argv[1].upper() # i.e. gkbs
channel = channel_dict[station[-1]]


# origin time (kaikoura earthquake)
origin = UTCDateTime('2016-11-13T11:02:56')
start = origin - 200
end = origin + 2000
year = start.year
julday = start.julday

# ===================== READ IN DATA =====================
mseed_files, resp_files = geonet_internal(station=station,
                                                channel=channel,
                                                start=origin,
                                                response=True)
st = Stream()
for files in mseed_files:
    st += read(files)
inv = read_inventory(resp_files[0])
for files in resp_files[1:]:
    inv += read_inventory(files)

# ===================== PREPROCESSING =====================
st.detrend('simple')
st.trim(starttime=start,endtime=end)
st.taper(max_percentage=0.05)
st.attach_response(inventories=inv)
st.remove_response(output='VEL',water_level=100)

# filter window
tmin = 1
tmax = 50
freqmin = 1/tmax
freqmax = 1/tmin
st.detrend("simple")
st.taper(max_percentage=0.05)
st.filter("bandpass",freqmin=freqmin,freqmax=freqmax,corners=3)


# ===================== SEPARATE DATA STREAMS =====================
north = st.select(component='1')
east = st.select(component='2')
if (len(north) and len(east)) == 0:
    north = st.select(component='N')
    east = st.select(component='E')
north = north[0].data
east = east[0].data

vertical = st.select(component='Z')[0].data
horizontal = np.sqrt(north**2 + east**2)
groundmotion = np.sqrt(vertical**2 + north**2 + east**2)
verticalsquared = np.sqrt(vertical**2)

stats = st[0].stats
samp_rate = stats.sampling_rate

# ===================== PLOTTING =====================
f,(ax1,ax1a,ax1b,ax3,ax3a,ax3b) = plt.subplots(6,
                                    sharex=True,
                                    sharey=False,
                                    figsize=(9,5))

# ===================== SUBPLOT 1,1a,1b (waveform) =====================
t = np.linspace(0,stats.endtime-stats.starttime,stats.npts)
ax1.plot(t,vertical,linewidth=1,c='k')
ax1.set_ylabel('Z (m/s)')
ax1a.plot(t,north,linewidth=1,c='k')
ax1a.set_ylabel('N/1 (m/s)')
ax1b.plot(t,east,linewidth=1,c='k')
ax1b.set_ylabel('E/2 (m/s)')
for ax in [ax1,ax1a,ax1b]:
    ax.grid(which="both")

plot_title = "Kaikoura Earthquake | "+"{instr} | {otime} | {t0}-{t1} s".format(
            otime=stats.starttime,
            instr=st[0].get_id()[:-4],
            t0=tmin,
            t1=tmax)
ax1.set_title(plot_title,fontsize=10)

# ===================== PROCESSING FOR SUBPLOT 3 =====================
duration,tover_plot,aover_plot,threshold_plot = [],[],[],[]
component_list = [verticalsquared,horizontal,groundmotion]

threshold_percentage = 0.1
for i,seismo in enumerate(component_list):
    # normalize trace between 0 and 1 for comparison between stations
    seismo /= seismo.max()
    seismo **= 2

    # find peak amplitude, set threshold
    peak_amp = seismo.max()
    threshold = peak_amp * threshold_percentage
    threshold_plot.append(threshold)

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

    duration.append(sum(sections))
    tover_plot.append(t_over)
    aover_plot.append(a_over)
    print("Duration criteria: {}".format(sum(sections)))

# ===================== SUBPLOT 3,3a,3b (duration criteria) ====================
axes = [ax3,ax3a,ax3b]
labels = ['Z','N+E','Z+N+E']
for AX,DU,CL,TO,AO,LA,TH in zip(
    axes,duration,component_list,tover_plot,aover_plot,labels,threshold_plot):

    # plot
    AX.plot(t,CL,'k')
    AX.scatter(TO,AO,c='r',marker='.',s=1,zorder=100)
    AX.set_ylabel('{}'.format(LA))
    # set threshold line and annotation
    h_lab = "Threshold = {}% peak amplitude".format(
                                                int(threshold_percentage*100))
    AX.axhline(y=TH,
                xmin=t[0],
                xmax=t[-1],
                zorder=1,
                color='r',
                linestyle='-.',
                linewidth=.85,
                label=h_lab)
    ano_x = TO[-1]
    ano_y = TH*2
    AX.annotate("Duration criteria: {}".format(round(DU,2)),
                                                xy=(ano_x,ano_y),
                                                xytext=(ano_x,ano_y))
    AX.grid()

ax3b.set_xlabel("Time (sec)")


# ===================== FINAL FIGURE ADJUSTMENTS =====================
plt.xlim([200,1500])
plt.subplots_adjust(wspace=.5, hspace=0)

figure_folder = '/seis/prj/fwi/bchow/spectral/output_plots/waveforms/kaikoura/'
figure_name = "{}.png".format(station)
outpath = os.path.join(figure_folder,figure_name)
f.savefig(outpath,dpi=250)
plt.show()
