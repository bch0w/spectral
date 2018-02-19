import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors

from random import shuffle
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
    colorrange = [scalarMap.to_rgba(i) for i in range(num_colors)]
    shuffle(colorrange)
    ax.set_prop_cycle('color',colorrange)

# naming conventions
# station_dict = {"RD01":"PRWZ", "RD02":"ANWZ", "RD03":"TURI", "RD04":"PORA",
#                 "RD05":"MNHR", "RD06":"DNVZ", "RD07":"WPAW", "RD08":"RAKW",
#                 "RD09":"MCHZ", "RD10":"CKHZ", "RD11":"KAHU", "RD12":"KWHZ",
#                 "RD13":"KERE", "RD14":"PNUI", "RD15":"WPUK", "RD16":"OROA",
#                 "RD17":"TEAC", "RD18":"RANC"}
# colocated short period sensors
station_dict = {"RD01":"PRWZ", "RD02":"ANWZ", "RD03":"", "RD04":"PRHZ",
                "RD05":"", "RD06":"DVHZ", "RD07":"", "RD08":"",
                "RD09":"MCHZ", "RD10":"CKHZ", "RD11":"KAHZ", "RD12":"KWHZ",
                "RD13":"KRHZ", "RD14":"PNHZ", "RD15":"WPHZ", "RD16":"",
                "RD17":"", "RD18":""}

month_dict = {"001*031":"JAN","032*059":"FEB","060*090":"MAR","091*120":"APR",
                "121*151":"MAY","152*181":"JUN","182*212":"JUL","213*243":"AUG",
                "244*273":"SEP","274*304":"OCT","305*334":"NOV","335-365":"DEC"}

# set path
npz_path = pathnames()['ppsd'] + "/RDF_decimateby5/" #temp array
# npz_files = glob.glob(npz_path + "*.npz")
rdf_files = glob.glob(npz_path + "RD*.npz")
ehz_files = glob.glob(npz_path + "*EHZ*.npz")

rdf_files.sort()
ehz_files.sort()

# start figure
f = plt.figure(dpi=200)
ax = f.add_subplot(111)

# unique colors and linestyles
cr = color_cycle(ax,len(rdf_files),'nipy_spectral')
line_styles = ['solid','dashed','dashdot','dotted']
num_styles = len(line_styles)
dashdotlist = [1,9,5,10,12,11,8,13,3,0,14]

for i,fid in enumerate(rdf_files):
    if i in dashdotlist:
        linestyle = 'dashed'
        linewidth = 0.5
        zorder = 100
    else:
        linestyle = 'dotted'
        linewidth = 0.5
        zorder = 2
    sta,cha,year = os.path.basename(fid).split(".")[:3]
    ppsd = PPSD.load_npz(fid)
    # avg = ppsd.get_mode()
    avg = ppsd.get_percentile(percentile=50)
    plt.plot(avg[0],avg[1],
                linestyle=linestyle,
                linewidth=linewidth,
                zorder=zorder,
                label="{}".format(sta))

colorrange = [cr[1],cr[9],cr[5],cr[10],cr[12],cr[11],cr[8],cr[13],cr[3],cr[0],cr[14]]
for i,fid in enumerate(ehz_files):
    sta,cha,year = os.path.basename(fid).split(".")[:3]
    ppsd = PPSD.load_npz(fid)
    # avg = ppsd.get_mode()
    avg = ppsd.get_percentile(percentile=50)

    plt.plot(avg[0],avg[1],
                linestyle='solid',
                linewidth=0.5,
                color=colorrange[i],
                zorder=3,
                label="{}".format(sta))

# plot lines for noise models and microseisms
nlnm_x,nlnm_y = get_nlnm()
nhnm_x,nhnm_y = get_nhnm()
plt.plot(nlnm_x,nlnm_y,'gray',alpha=0.7,zorder=1)
plt.plot(nhnm_x,nhnm_y,'gray',alpha=0.7,zorder=1)
# plt.axvline(x=16,color='k',alpha=0.5,zorder=1,lw=0.5,ls='dashed')
# plt.axvline(x=4,color='k',alpha=0.5,zorder=1,lw=0.5,ls='dashed')
# plt.axvline(x=30,color='k',alpha=0.5,zorder=1,lw=0.5,ls='dashed')
# ax.axvspan(3, 10, alpha=0.1, color='green',zorder=1)

plt.legend(ncol=3,prop={"size":5})

# for geonet permanent stations
plt.ylim([nlnm_y.min(),-90])
plt.xlim([0.2,100])

# for rdf temp stations
# plt.xlim([0.2,100])
# plt.ylim([nlnm_y.min(),-90])

plt.xscale("log")
plt.xlabel("Period (s)")
plt.ylabel("Amplitude [m^2/s^4/Hz][dB]")
# plt.title("Mean values of year-long PPSD\'s for GEONET permanent seismometers\n"
#             "Year: 2015 | Sampling Rate: 10 Hz | # Stations: {}".format(len(npz_files)))
plt.title("Median values of PPSD for RDF stations and colocated short-period sensors")
plt.grid()
# plt.savefig(figure_savename)
plt.show()

# plot variations from baseline
    # f2 = plt.figure(dpi=200)
    # # plt.plot(avg[0],noncoastal_mean,label="Noncoastal mean",linewidth=0.75)
    # for coast,sta in zip(coastal_avgs,coastal_sta):
    #     difference = np.subtract(coast,noncoastal_mean)
    #     plt.plot(avg[0],difference,label=sta,linewidth=0.75)
    # for noncoast in noncoastal_avgs:
    #     difference = np.subtract(noncoast,noncoastal_mean)
    #     plt.plot(avg[0],difference,linewidth=0.25,alpha=0.25,color="gray",zorder=0)
    #
    # plt.xscale("log")
    # plt.xlabel("Period (s)")
    # plt.ylabel("Ampltidue [dB]")
    # # yticks = range(-12,13,2)
    # # plt.yticks(yticks)
    # # plt.ylim([-12,12])
    # # plt.title("Variation of non-coastal average and coastal site PPSD")
    # plt.title("Variation of RDF array quiet vs noisy station PPSD averages")
    # plt.grid()
    # plt.legend(loc='best',prop={'fontsize':0.05})
    # plt.show()
