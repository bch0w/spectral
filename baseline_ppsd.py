"""create baseline variation plots for basin vs nonbasin stations
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
    # shuffle(colorrange)
    ax.set_prop_cycle('color',colorrange)

#  =================================== MAIN ===================================
station_list = ['PUZ', 'TSZ', 'TLZ', 'HAZ', 'KNZ', 'TOZ', 'KHEZ', 'MKAZ', 'MRZ',
                'BFZ', 'WAZ', 'RTZ', 'KUZ', 'GRZ', 'BKZ', 'WCZ', 'URZ',
                'MXZ', 'PXZ', 'OUZ', 'MWZ', 'RATZ', 'OPRZ', 'HIZ']
                #'VRZ', 'WSRZ'

npz_path = os.path.join(pathnames()['ppsd'],"geonet_db5_averages","")
periods = np.load(npz_path + "periods.npy")

for component in ["Z","N","E","H"]:
    npz_files = glob.glob(npz_path + "*HH{}*".format(component))

# =========================== PLOT 1 BASIN/NONBASIN ===========================
    f,(ax1,ax2) = plt.subplots(2,sharex=True,sharey=False,
                                dpi=200,figsize=(7,7))

    basin = ["KNZ","PUZ"]
    coastal = ["BFZ","PXZ","KNZ","PUZ","MXZ"] #actually coastal
    color_dict = {"BFZ":"r",
                "PXZ":"deeppink",
                "KNZ":"purple",
                "PUZ":"orange",
                "MXZ":"c"}
    basin_sta,basin_avg,nonbasin_avg = [],[],[]
    # color_cycle(ax1,len(coastal),'plasma')

    for i,fid in enumerate(npz_files):
        sta = os.path.basename(fid).split('.')[0]
        # skip over noisy stations
        if sta not in station_list:
            print(sta)
            continue
        avg = np.load(fid)
        if sta in coastal:
            basin_avg.append(avg)
            basin_sta.append(sta)
            if sta in basin:
                linestyle,alpha,zorder = 'solid',1,3
            else:
                linestyle,alpha,zorder = 'dashed',1,3

            ax1.plot(periods,avg,
                        linestyle=linestyle,
                        linewidth=1,
                        color=color_dict[sta],
                        alpha=alpha,
                        zorder=zorder,
                        label="{}".format(sta))
        else:
            nonbasin_avg.append(avg)
            ax1.plot(periods,avg,
                        linestyle='dashed',
                        linewidth=0.6,
                        alpha=0.2,
                        color='b',
                        zorder=1)


    # take averages of all lines and plot
    average_of_basin_avg = np.median(np.array(basin_avg),axis=0)
    average_of_nonbasin_avg = np.median(np.array(nonbasin_avg),axis=0)
        # average_of_basin_avg = np.array(basin_avg).mean(axis=0)
        # average_of_nonbasin_avg = np.array(nonbasin_avg).mean(axis=0)
    ax1.plot(periods,average_of_nonbasin_avg,
                    linestyle='solid',
                    linewidth=1.55,
                    alpha=1,
                    color='b',
                    zorder=4,
                    label="non-coastal median")

    # plot low and high noise models
    nlnm_x,nlnm_y = get_nlnm()
    nhnm_x,nhnm_y = get_nhnm()
    ax1.plot(nlnm_x,nlnm_y,'k',alpha=.8,zorder=1,linewidth=0.35)
    ax1.plot(nhnm_x,nhnm_y,'k',alpha=.8,zorder=1,linewidth=0.35)

    # plot parameters
    ax1.set_title("2015 Geonet Permanent Seismometers\n"
                "Coastal vs. Noncoastal Stations {}".format(component))
    ax1.set_ylabel("Amplitude [m^2/s^4/Hz][dB]")

    ax1.set_ylim([-165,-96])
    ax1.set_xlim([1,100])
    ax1.set_xscale("log")

    ax1.set_axisbelow(True)
    ax1.tick_params(which='both',direction='in',top=True,right=True)
    ax1.minorticks_on()
    ax1.grid(which='minor',linestyle=':',linewidth='0.5',color='k',alpha=0.25)
    ax1.grid(which='major',linestyle='-',linewidth='0.5',color='k',alpha=0.15)

    ax1.legend(ncol=2,prop={"size":5})

# =========================== PLOT 2 VARIATIONS ==========================
    # plot variations from baseline
    color_cycle(ax2,len(coastal),'plasma')
    for coast,sta in zip(basin_avg,basin_sta):
        difference = np.subtract(coast,average_of_nonbasin_avg)
        if sta in basin:
            linestyle = 'solid'
        else:
            linestyle = 'dashed'
        ax2.plot(periods,difference,
                    color=color_dict[sta],
                    label=sta,
                    linewidth=0.75,
                    linestyle=linestyle)
    for noncoast in nonbasin_avg:
        difference = np.subtract(noncoast,average_of_nonbasin_avg)
        ax2.plot(periods,difference,
                    linewidth=0.25,
                    alpha=0.25,
                    color="gray",
                    zorder=1)

    # doesnt work?
    ax2.hlines(y=0,xmin=1,xmax=100,
                alpha=.9,linestyle='solid',color='k',linewidth=0.5,zorder=1)

    ax2.set_xscale("log")
    ax2.set_xlabel("Period (s)")
    ax2.set_ylabel("Ampltidue [dB]")

    ax2.set_axisbelow(True)
    ax2.tick_params(which='both',direction='in',top=True,right=True)
    ax2.minorticks_on()
    ax2.grid(which='major',linestyle='-',linewidth='0.5',color='k',alpha=0.15)
    ax2.grid(which='minor',linestyle=':',linewidth='0.5',color='k',alpha=0.25)

    ax2.set_ylim([-10,15])
    # ax2.set_title("Variation of Coastal and Noncoastal baseline")
    ax2.set_xlim([1,100])
    ax2.legend(ncol=2,prop={"size":5})


    # final plot adjustments
    plt.subplots_adjust(wspace=.5, hspace=.05)
    plt.tight_layout
    plt.savefig(pathnames()['plots']+
                        "ppsd_plots/variations{}.png".format(component))
    plt.show()
