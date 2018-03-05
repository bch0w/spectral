"""19.2.18
PARED AND EDITED SCRIPT TO MAKE EQ DURATION PLOTS FOR VERTICAL COMPONENT
FIGURES FOR SLOW SLIP WORKSHOP PRESENTATION

Plot waveforms of all three components for a given GEONET permanent station,
or a temporary RDF station, with preprocessing and filtering set in the script.
Also subplots of determining duration criteria, which is captured using an
amplitude threshold criteria and summing up the time sections (dt) where the
waveform crosses this amplitude threshold.

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
from obspy.clients.fdsn import Client
from getdata import pathnames, event_stream

# ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# global plot parameters
# mpl.rcParams['text.usetex'] = True
# mpl.rcParams['font',**{'family':'sans-serif','sans-serif':['Helvetica']}]
mpl.rcParams['font.size'] = 15
mpl.rcParams['lines.linewidth'] = 1


# ====================================== MAIN ==================================
# set station and channel ID
channel_dict = {"Z":"HH*","S":"BN*"}
station = sys.argv[1].upper() # i.e. gkbs
event_id = sys.argv[2] #2014p240655,2015p822263,2016p892721,2017p059122
channel = channel_dict[station[-1]]

# for determining amplitude threshold
threshold_percentage = 0.125

# origin = UTCDateTime('2016-11-13T11:02:56') # kaikoura

# ================================ READ IN DATA ===============================
pad = 200
st,inv,cat = event_stream(station=station,
                            channel="HHZ",
                            event_id=event_id,
                            pad=pad)
# set timing
event = cat[0]
origin = event.origins[0].time
start = origin - pad
end = origin + 2000
year = start.year
julday = start.julday

# ================================ PREPROCESSING ===============================
st.detrend('simple')
st.trim(starttime=start,endtime=end)
st.detrend('linear')
st.taper(max_percentage=0.05)
st.attach_response(inventories=inv)
st.remove_response(output='VEL',water_level=100)#,plot=True)
plt.show()

# filter window
tmin = 5
tmax = 30
freqmin = 1/tmax
freqmax = 1/tmin
st.detrend("simple")
st.taper(max_percentage=0.05)
st.filter("bandpass",freqmin=freqmin,freqmax=freqmax,corners=3)

tracestart = origin
traceend = origin + 700
st.trim(starttime=tracestart,endtime=traceend)

# ================================ SEPARATE DATA STREAMS =======================
vertical = st.select(component='Z')[0].data
vertical *= (10**3)

stats = st[0].stats
samp_rate = int(stats.sampling_rate)
t = np.linspace(0,stats.endtime-stats.starttime,stats.npts)


# ================================ PLOTTING ====================================
f,(ax1,ax3) = plt.subplots(2,
                            sharex=True,
                            sharey=False,
                            figsize=(9,5))

# ================================ SUBPLOT 1,1a,1b (waveform) ==================
ax1.plot(t,vertical,linewidth=1.5,c='k')
ax1.set_ylabel('Vertical Velocity (mm/s)')
ax1.grid(which="both")

# ===================== AMPLITUDE PROCESSING FOR SUBPLOT 3 =====================
duration_a,tover_plot,aover_plot,threshold_plot,sample_plot = [],[],[],[],[]

# for each component, vertical, horizontal and vector sum
seismo=np.sqrt(vertical**2)
peak_amp = seismo.max()
threshold = peak_amp * threshold_percentage
threshold_plot.append(threshold)

# loop over seismogram, determine start and end of peak energy
# a for amplitude, s for sample
a_over, s_over = [],[]
for i,amplitude in enumerate(seismo):
    if amplitude >= threshold:
        s_over.append(i)
        a_over.append(amplitude)

# find edgepoints by checking if the next sample j is the same as i+1
s_edge,a_edge,sections,samples = [s_over[0]],[a_over[0]],[],[]
for i,(S,A) in enumerate(zip(s_over[1:-2],a_over[1:-2])):
    if s_over[i+2] != (S + 1):
        section = np.trapz(a_edge,s_edge)
        sections.append(section)

        # determine number of samples covered
        samples.append(len(s_edge))
        s_edge,a_edge = [s_over[i+1]],[a_over[i+1]]
    else:
        # if the next sample is the same, keep going
        s_edge.append(S)
        a_edge.append(A)

# convert samples to time
t_over = []
for S in s_over:
    t_over.append(t[S])

# create lists for plotting
sample_plot.append(sum(samples)/samp_rate)
duration_a.append(sum(sections))
tover_plot.append(t_over)
aover_plot.append(a_over)

# ==================== SUBPLOT 3,3a,3b (amplitude criteria) ====================
axes = [ax3]
component_list = [seismo]

labels = ['']
for AX,DU,CL,TO,AO,LA,TH,SA in zip(axes,duration_a,component_list,
                                tover_plot,aover_plot,labels,
                                threshold_plot,sample_plot):

    # plot
    AX.plot(t,CL,'k')
    AX.scatter(TO,AO,c='r',marker='x',s=0.2,zorder=100)
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
                linewidth=1.5,
                label=h_lab)
    ano_x = TO[-1]
    ano_y = TH*2
    # AX.annotate("Duration criteria: {}s".format(round(SA,2)),
    #                                             xy=(ano_x,ano_y),
    #                                             xytext=(ano_x,ano_y))
    AX.grid()

ax3.set_xlabel("Time (sec)")

# ====================== FINAL FIGURE ADJUSTMENTS ==============================
plot_title = "{eid} | {instr} | {otime} | {t0}-{t1}s | {T} | {D}s | {P}".format(
            eid=event_id,
            otime=stats.starttime,
            instr=st[0].get_id()[:-4],
            t0=tmin,
            t1=tmax,
            T="{}%".format(round(threshold_percentage * 100,2)),
            D=round(SA,2),
            P=peak_amp)

ax1.set_title(plot_title,fontsize=10)
plt.xlim([0,500])
plt.subplots_adjust(wspace=.5, hspace=0)

figure_folder = pathnames()["plots"]+ 'waveforms/{}/'.format(event_id)
if not os.path.exists(figure_folder):
    os.makedirs(figure_folder)
figure_name = "{0}_{1}-{2}.png".format(station,tmin,tmax)
outpath = os.path.join(figure_folder,figure_name)


# f.savefig(outpath,dpi=250)
# plt.show()

with open(figure_folder + '{}_{}-{}_amplitudes.txt'.format(event_id,tmin,tmax), 'a+') as f:
    f.write('{0} {1}\n'.format(station,peak_amp))
