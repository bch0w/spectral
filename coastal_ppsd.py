"""create baseline variation plots for coastal vs noncoastal stations
for the geonet permanent broadband seismometer network
"""
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors

from random import shuffle
from getdata import pathnames, vog
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

vic_or_gns = vog()

# set path
npz_path = os.path.join(pathnames(vic_or_gns)['ppsd'],
                        "geonet_db5_averages","")
station_list = ['PUZ', 'TSZ', 'TLZ', 'HAZ', 'KNZ', 'TOZ', 'KHEZ', 'MKAZ', 'MRZ',
                'BFZ', 'WAZ', 'RTZ', 'WCZ', 'WSRZ', 'KUZ', 'GRZ', 'BKZ',
                'MXZ', 'PXZ', 'OUZ', 'MWZ', 'RATZ', 'OPRZ', 'HIZ'] #VRZ, URZ

choice = 'Z'
npz_files = glob.glob(npz_path + "*HH{}*".format(choice))

# start figure
f = plt.figure(dpi=200)
ax = f.add_subplot(111)

# unique colors and linestyles
cr = color_cycle(ax,len(station_list),'nipy_spectral')
linestyles = ['solid','dashed','dashdot','dotted']
coastal = ["BFZ","PXZ","KNZ","PUZ","MXZ"]

periods = np.load(npz_path + "periods.npy")
for i,fid in enumerate(npz_files):
    sta = os.path.basename(fid).split('.')[0]
    if sta in coastal:
        c='b'
    else:
        c='k'
    avg = np.load(fid)
    plt.plot(periods,avg,
                linestyle=linestyles[i%len(linestyles)],
                linewidth=0.5,
                # color=c,
                zorder=3,
                label="{}".format(sta))

# plot lines for noise models and microseisms
nlnm_x,nlnm_y = get_nlnm()
nhnm_x,nhnm_y = get_nhnm()
plt.plot(nlnm_x,nlnm_y,'gray',alpha=0.7,zorder=1)
plt.plot(nhnm_x,nhnm_y,'gray',alpha=0.7,zorder=1)

plt.legend(ncol=3,prop={"size":5})

# for geonet permanent stations
plt.ylim([nlnm_y.min(),-90])
plt.xlim([0.2,100])

plt.xscale("log")
plt.xlabel("Period (s)")
plt.ylabel("Amplitude [m^2/s^4/Hz][dB]")

plt.title("PPSD Averages for Geonet permanent seismometers | {}".format(choice))
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





# average monthly data by station and component
    # npz_path = os.path.join(pathnames(vic_or_gns)['ppsd'],
    #                         "geonet_monthly_decimateby5",
    #                         "2015",
    #                         "*",
    #                         "")
    # station_list = ['PUZ', 'TSZ', 'TLZ', 'HAZ', 'KNZ', 'TOZ', 'KHEZ', 'MKAZ', 'MRZ',
    #                 'BFZ', 'WAZ', 'RTZ', 'WCZ', 'WSRZ', 'KUZ', 'GRZ', 'BKZ',
    #                 'MXZ', 'PXZ', 'OUZ', 'MWZ', 'RATZ', 'OPRZ', 'HIZ'] #VRZ, URZ
    # all_averages = []
    # for station in station_list:
    #     station_averages = []
    #     for comp in ['N','E','Z']:
    #         list_of_files = glob.glob(npz_path + "{s}*HH{c}*.npz".format(s=station,
    #                                                                   c=comp))
    #         # print(station,comp,len(list_of_files))
    #         average_avg = []
    #         for fname in list_of_files:
    #             ppsd_temp = PPSD.load_npz(fname)
    #             avg_temp = ppsd_temp.get_mode()
    #             average_avg.append(avg_temp[1])
    #
    #         mean_of_averages = np.array(average_avg).mean(axis=0)
    #         station_averages.append(mean_of_averages) #should be 3x25
    #
    #     all_averages.append(station_averages)
    #
    # periods = avg_temp[0]
    #
    # for i,AVG in enumerate(all_averages):
    #     north,east,vert = AVG
    #     hori = np.array([north,east]).mean(axis=0)
    #     np.save("./ppsd_arrays/geonet_db5_averages/{}.HHZ.2015".format(station_list[i]),vert)
    #     np.save("./ppsd_arrays/geonet_db5_averages/{}.HHN.2015".format(station_list[i]),north)
    #     np.save("./ppsd_arrays/geonet_db5_averages/{}.HHE.2015".format(station_list[i]),east)
    #     np.save("./ppsd_arrays/geonet_db5_averages/{}.HHH.2015".format(station_list[i]),hori)
