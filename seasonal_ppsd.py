"""Create PPSD's for geonet permanent stations by seasonal binning
"""
import os
import sys
sys.path.append('./..')
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm

from getdata import pathnames
from obspy.signal import PPSD
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm

def color_cycle(ax,length,cmap):
    """sets up a cycle of colors to be used in figures with many lines
    """
    import matplotlib.colors as colors
    num_colors = length
    cm = plt.get_cmap(cmap)
    norm_col = colors.Normalize(vmin=0,vmax=num_colors-1)
    scalarMap = mplcm.ScalarMappable(norm=norm_col,cmap=cm)
    colors = [scalarMap.to_rgba(i) for i in range(num_colors)]
    ax.set_prop_cycle('color',colors)

# main
vic_or_gns = "vic"

month_dict = {"JAN":"001*031",
              "FEB":"032*059",
              "MAR":"060*090",
              "APR":"091*120",
              "MAY":"121*151",
              "JUN":"152*181",
              "JUL":"182*212",
              "AUG":"213*243",
              "SEP":"244*273",
              "OCT":"274*304",
              "NOV":"305*334",
              "DEC":"335*365"}

station_list = ["TSZ","MRZ","BFZ","WAZ","KHEZ","VRZ","HIZ","TLZ","TOZ",
                "KUZ","GRZ","WSRZ","RATZ","WCZ","OUZ","MWZ","BKZ","KNZ","RTZ",
                "HAZ","PUZ","OPRZ","URZ","PXZ","MXZ"] #"MKAZ" < no files

npz_path = pathnames(vic_or_gns)['ppsd']

# average summer/winter months into separate PPSDs
winter_avgs,summer_avgs = [],[]
for station in station_list:
    winter_files,summer_files = [],[]
    # =========================== WINTER ===========================
    for winter in ["JUN","JUL","AUG"]:
        winter_files += glob.glob(npz_path + "{s}*{w}.npz".format(s=station,
                                                        w=month_dict[winter]))
    # grab mean/mode from winter
    average_avg = []
    for temp in winter_files:
        ppsd_temp = PPSD.load_npz(temp)
        avg_temp = ppsd_temp.get_mode()
        average_avg.append(avg_temp[1])

    # take mean of all winter months
    mean_of_averages = np.array(average_avg).mean(axis=0)
    winter_avgs.append(mean_of_averages)

    # =========================== SUMMER ===========================
    for summer in ["DEC","JAN","FEB"]:
        summer_files += glob.glob(npz_path + "{s}*{w}.npz".format(s=station,
                                                        w=month_dict[summer]))
    # grab mean/mode from each season
    average_avg = []
    for temp in summer_files:
        ppsd_temp = PPSD.load_npz(temp)
        avg_temp = ppsd_temp.get_mode()
        average_avg.append(avg_temp[1])

    # take mean of all summer months
    mean_of_averages = np.array(average_avg).mean(axis=0)
    summer_avgs.append(mean_of_averages)

import ipdb;ipdb.set_trace()

# PPSD x-axis; same for all
periods = avg_temp[0]

# find averages for winter and summer
average_winter = np.array(winter_avgs).mean(axis=0)
average_summer = np.array(summer_avgs).mean(axis=0)

# start figure
f = plt.figure(dpi=200)
ax = f.add_subplot(111)
line_styles = ['solid','dashed','dashdot','dotted']
num_styles = len(line_styles)

color_cycle(ax,len(winter_avgs),'winter')
for i,(winter,sta) in enumerate(zip(winter_avgs,station_list)):
    plt.plot(periods,winter,linewidth=0.5,
                            alpha=0.5,
                            linestyle=line_styles[i%num_styles],
                            label=sta)

color_cycle(ax,len(summer_avgs),'autumn')
for i,(summer,sta) in enumerate(zip(summer_avgs,station_list)):
    plt.plot(periods,summer,linewidth=0.5,
                            alpha=0.5,
                            linestyle=line_styles[i%num_styles],
                            label=sta)

# plot average
plt.plot(periods,average_winter,color='c',
                                linestyle='solid',
                                label='WINT')
plt.plot(periods,average_summer,color='r',
                                linestyle='solid',
                                label='SUMM')



nlnm_x,nlnm_y = get_nlnm()
nhnm_x,nhnm_y = get_nhnm()
plt.plot(nlnm_x,nlnm_y,'gray',alpha=0.7)
plt.plot(nhnm_x,nhnm_y,'gray',alpha=0.7)

# microseism references
plt.axvline(x=16,color='k',alpha=0.5,zorder=1,lw=0.5,ls='dashed')
plt.axvline(x=4,color='k',alpha=0.5,zorder=1,lw=0.5,ls='dashed')
plt.axvline(x=30,color='k',alpha=0.5,zorder=1,lw=0.5,ls='dashed')

handles, labels = ax.get_legend_handles_labels()
# sort both labels and handles by labels
labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
ax.legend(handles, labels,ncol = 4,prop={"size":5})

# plt.legend(ncol=4,prop={"size":5})
plt.ylim([nlnm_y.min(),-90])
plt.xlim([0.2,100])
plt.xscale("log")
plt.xlabel("Period (s)")
plt.ylabel("Amplitude [m^2/s^4/Hz][dB]")
plt.title("Seasonal averages of PPSD\'s for GEONET permanent seismometers\n"
            "Year: 2015 | # Stations: {}".format(len(station_list)))
# plt.title("Mode values of PPSD for RDF Temporary Array")
plt.grid()
plt.show()
