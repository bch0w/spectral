import os
import sys
sys.path.append('./..')
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors

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

vic_or_gns = "vic"

station_dict = {"RD01":"PRWZ", "RD02":"ANWZ", "RD03":"TURI", "RD04":"PORA",
                "RD05":"MNHR", "RD06":"DNVZ", "RD07":"WPAW", "RD08":"RAKW",
                "RD09":"MCHZ", "RD10":"CKHZ", "RD11":"KAHU", "RD12":"KWHZ",
                "RD13":"KERE", "RD14":"PNUI", "RD15":"WPUK", "RD16":"OROA",
                "RD17":"TEAC", "RD18":"RANC"}

month_dict = {"001*031":"JAN","032*059":"FEB","060*090":"MAR","091*120":"APR",
                "121*151":"MAY","152*181":"JUN","182*212":"JUL","213*243":"AUG",
                "244*273":"SEP","274*304":"OCT","305*334":"NOV","335-365":"DEC"}

# set path
npz_path = pathnames(vic_or_gns)['ppsd'] + "/RDF_decimateby5/" #temp array
npz_files = glob.glob(npz_path + "RD*.npz")
# npz_files = glob.glob(npz_path + "*HHZ*001-365.npz")

# manual set npz file load to compare two stations (kidnapper)
npz_files = []
npz_files.append("/seis/prj/fwi/bchow/spectral/ppsd_arrays/RDF_decimateby5/RD10.HHZ.2017.321-2017.321.npz")
npz_files.append("/seis/prj/fwi/bchow/spectral/ppsd_arrays/CKHZ.EHZ.2017.321-2018.014.npz")
npz_files.append("/seis/prj/fwi/bchow/spectral/ppsd_arrays/RDF_decimateby5/RD11.HHZ.2017.321-2017.314.npz")
npz_files.append("/seis/prj/fwi/bchow/spectral/ppsd_arrays/KAHZ.EHZ.2017.321-2018.014.npz")


npz_files.sort()

# start figure
f = plt.figure(dpi=200)
ax = f.add_subplot(111)

# unique colors and linestyles
# color_cycle(ax,len(npz_files),'nipy_spectral')
# line_styles = ['solid','dashed','dashdot','dotted']
line_styles = ['solid','solid','dashed','dashed']
num_styles = len(line_styles)

color_styles = ['r','k','r','k']

# coastal_stations = ["BFZ", "PXZ", "KNZ", "PUZ"]
coastal_stations = ["RD12","RD03","RD11","RD10","RD17","RD16"]
coastal_avgs,coastal_sta,noncoastal_avgs = [],[],[]

for i,fid in enumerate(npz_files):
    sta,cha,year = os.path.basename(fid).split(".")[:3]
    ppsd = PPSD.load_npz(fid)
    # avg = ppsd.get_mode()
    avg = ppsd.get_percentile(percentile=50)
    if sta == "WSRZ":
        continue
    plt.plot(avg[0],avg[1],
                linestyle=line_styles[i%num_styles],
                color=color_styles[i%num_styles],
                linewidth=0.75,
                label="{}".format(sta))
    # if sta not in coastal_stations:
    #     noncoastal_avgs.append(avg[1])
    #     plt.plot(avg[0],avg[1],
    #             color='gray',
    #             alpha=0.25,
    #             zorder=0,
    #             linewidth=0.25)
    # else:
    #     coastal_sta.append(sta)
    #     coastal_avgs.append(avg[1])
    #     plt.plot(avg[0],avg[1],
    #             label="{sta}".format(sta=sta),
    #             linestyle=line_styles[i%num_styles],
    #             linewidth=.75)


# TEMPORARY
ppsd = PPSD.load_npz("/seis/prj/fwi/bchow/spectral/ppsd_arrays/RDF_decimateby5/RD06.HHZ.2017.320-2017.320.npz")
avg = ppsd.get_percentile(percentile=50)
sta= "RD06"
plt.plot(avg[0],avg[1],
            linestyle="dashed",
            color="orange",
            linewidth=0.75,
            label="{}".format(sta))
# TEMPORARY
# noncoastal_mean = np.array(noncoastal_avgs).mean(axis=0)
# plt.plot(avg[0],noncoastal_mean,linewidth=1,label="Quiet station average")

# plot lines for noise models and microseisms
nlnm_x,nlnm_y = get_nlnm()
nhnm_x,nhnm_y = get_nhnm()
plt.plot(nlnm_x,nlnm_y,'gray',alpha=0.7)
plt.plot(nhnm_x,nhnm_y,'gray',alpha=0.7)
plt.axvline(x=16,color='k',alpha=0.5,zorder=1,lw=0.5,ls='dashed')
plt.axvline(x=4,color='k',alpha=0.5,zorder=1,lw=0.5,ls='dashed')
plt.axvline(x=30,color='k',alpha=0.5,zorder=1,lw=0.5,ls='dashed')

plt.legend(ncol=1,prop={"size":5})

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
plt.title("Median values of PPSD for RDF stations vs colocated geonet short period sensors")
plt.grid()
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
