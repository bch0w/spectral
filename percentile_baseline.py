"""same as baseline_ppsd, but instead of looking only at coastal vs non coastal
stations, looks at all stations that lie about a certain percentile threshold
"""
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors
from obspy.signal import PPSD
from random import shuffle
from getdata import pathnames
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
station_list = ['PUZ', 'TSZ', 'TLZ', 'HAZ', 'KNZ', 'TOZ', 'KHEZ', 'MRZ',
                'BFZ', 'WAZ', 'RTZ', 'KUZ', 'GRZ', 'BKZ', 'WCZ', 'URZ',
                'MXZ', 'PXZ', 'OUZ', 'MWZ', 'RATZ', 'OPRZ', 'HIZ']
                #'VRZ', 'WSRZ' 'MKAZ'

npz_path = os.path.join(pathnames()['ppsd'],"geonet_db5_summer1718","")
periods = np.load(npz_path + "periods.npy")

for component in ["Z","N","E","H"]:
    npz_files = glob.glob(npz_path + "*HH{}*".format(component))

# ===================== PRCOESSING 1 PERCENTILE ==============================
    # put all data into two lists
    station_list, ppsd_list = [],[]
    for i,fid in enumerate(npz_files):
        sta = os.path.basename(fid).split('.')[0]
        station_list.append(sta)
        # skip over noisy stations
        if sta not in station_list:
            print(sta)
            continue
        _,ext = os.path.splitext(fid)
        # read in data, could be npy or npz files
        if ext == ".npy":
            avg = np.load(fid)
        elif ext == ".npz":
            ppsd = PPSD.load_npz(fid)
            avg = ppsd.get_percentile(percentile=50)[1]
        ppsd_list.append(avg)

    # find percentile to compare all other stations
    p25 = np.percentile(a=np.array(ppsd_list),q=25,axis=0)
    p50 = np.percentile(a=np.array(ppsd_list),q=50,axis=0)
    p75 = np.percentile(a=np.array(ppsd_list),q=75,axis=0)
    percentiles = [p25,p50,p75]

# =========================== PLOT 1 BASIN/NONBASIN ===========================

    # f,(ax1,ax2) = plt.subplots(2,sharex=True,sharey=False,
    #                             dpi=200,figsize=(7,7))
    f,ax1 = plt.subplots(1,dpi=200,figsize=(3,3))
    color_cycle(ax1,len(ppsd_list),'nipy_spectral')
    for avg,sta in zip(ppsd_list,station_list):
        ax1.plot(periods,avg,
                    linestyle="solid",
                    linewidth=1,
                    alpha=1,
                    zorder=1,
                    label="{}".format(sta))

    for p in percentiles:
        ax1.plot(periods,p,
                    linestyle="dashed",
                    linewidth=2,
                    color='k',
                    zorder=2,
                    label='75%')

    # plot low and high noise models
    nlnm_x,nlnm_y = get_nlnm()
    nhnm_x,nhnm_y = get_nhnm()
    ax1.plot(nlnm_x,nlnm_y,'k',alpha=.8,zorder=1,linewidth=0.35)
    ax1.plot(nhnm_x,nhnm_y,'k',alpha=.8,zorder=1,linewidth=0.35)

    # plot parameters
    ax1.set_title("Summer 17/18 Geonet Permanent Seismometers\n"
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
    plt.show()
    sys.exit()

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

    ax2.hlines(y=0,xmin=1,xmax=100,
                alpha=.9,linestyle='solid',color='k',linewidth=0.5,zorder=1)

    ax2.set_xscale("log")
    ax2.set_xlabel("Period (s)")
    ax2.set_ylabel("Baseline Variation [dB]")

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
                    "ppsd_plots/summer1718_variations{}.png".format(component))
    # plt.show()
