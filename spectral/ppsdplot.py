import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal import PPSD
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm

sys.path.append('../modules')
from getdata import pathnames
from plotmod import pretty_grids, build_color_dictionary

import matplotlib as mpl
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1.5

# ============================== HELPER FUNCTIONS ==============================
def station_dict(code):
    station_dict = {"RD01":"PRWZ", "RD02":"ANWZ", "RD03":"TURI", "RD04":"PORA",
                    "RD05":"MNHR", "RD06":"DNVZ", "RD07":"WPAW", "RD08":"RAKW",
                    "RD09":"MCNL", "RD10":"CKHZ", "RD11":"KAHU", "RD12":"KWHZ",
                    "RD13":"KERE", "RD14":"PNUI", "RD15":"WPUK", "RD16":"OROA",
                    "RD17":"TEAC", "RD18":"RANC", "RD19":"MATT","RD20":"KAHU2",
                    "RD21":"TEAC2"}
    return station_dict[code]

def prime_plot():
    """set up plot objects
    """
    f = plt.figure(figsize=(9,5),dpi=200)
    ax = f.add_subplot(111)
    pretty_grids(ax)

    # plot lines for noise models and microseisms
    nlnm_x,nlnm_y = get_nlnm()
    nhnm_x,nhnm_y = get_nhnm()
    plt.plot(nhnm_x,nhnm_y,'gray',alpha=0.1,linewidth=1)
    plt.plot(nlnm_x,nlnm_y,'gray',alpha=0.1,linewidth=1)
    ax.fill_between(nhnm_x,nlnm_y,nhnm_y,facecolor='gray',alpha=0.1)

    # set common plotting parameters
    plt.xlim([0.2,100])
    plt.ylim([nlnm_y.min(),-90])
    plt.xscale("log")
    plt.xlabel("Period [s]")
    plt.ylabel("Amplitude [m^2/s^4/Hz][dB]")

    return f, ax

# ========================== PROCESSING FUNCTIONS ==============================
def multi_component_ppsd(comp,specific,avg='median'):
    """Return average arrays for each station in a given folder, choice of
    component, can be N/E/Z/*
    Choice of average available from the PPSD objects
    """
    fidtmplt = "RD{s:0>2}.HH{c}*"

    if specific == "ALL":
        folder = 'FATHOM/*'
    else:
        folder = 'FATHOM/{}'.format(specific)
    folderlist = glob.glob(os.path.join(pathnames()['ppsd'],folder))

    averages,stations = [],[]
    for F in folderlist:
        print(F)
        # check range of stations available using filenames
        allfiles = glob.glob(os.path.join(F,'*HH*.npz'))
        allfiles.sort()
        lowrange = int(os.path.basename(allfiles[0]).split('.')[0][2:])
        highrange = int(os.path.basename(allfiles[-1]).split('.')[0][2:])

        # loop by station and inner loop by components per station
        for i in range(lowrange,highrange+1):
            files = glob.glob(os.path.join(F,fidtmplt.format(s=i,c=comp)))
            if not files:
                print('no files',fidtmplt.format(s=i,c=comp))
                continue
            inneraverages = []
            for fid in files:
                ppsd = PPSD.load_npz(fid)
                if avg == 'mean':
                    average = ppsd.get_mean()
                elif avg == 'mode':
                    average = ppsd.get_mode()
                elif avg  == 'median':
                    average = ppsd.get_percentile(percentile=50)
                inneraverages.append(average[1])

            outeraverage = np.mean(np.array(inneraverages),axis=0)
            averages.append(outeraverage)
            stations.append("RD{:0>2}".format(i))
    periods = average[0]

    if specific == "ALL":
        stations,averages = collapse_stations(stations,averages)

    return stations,periods,averages

def collapse_stations(stations,averages):
    """if stations are repeated (e.g. from muilticomponent ppsd gather) take
    the mean of all arrays and return collapsed arrays
    """
    stations = np.array(stations)
    averages = np.array(averages)
    setstations = list(set(stations))
    setstations.sort()
    stationsout,averagesout = [],[]
    for sta in setstations:
        indices = np.where(stations==sta)[0]
        avginner = []
        for i in indices:
            avginner.append(averages[i])
        avgouter = np.mean(avginner,axis=0)
        averagesout.append(avgouter)
        stationsout.append(sta)

    return stationsout,averagesout


# ============================== PLOTTING FUNCTIONS ============================
def plot_ppsd(save=False,show=True):
    """Main plotting function, plots average lines from a given folder with
    various distinct colors and markers.

    ++ parameters can be the following:
    specific: 'ALL','JULY17_JAN18','JAN18_MAR18','MAR18_MAY18','MAY18_JUNE18'
    component: 'N','E','Z','*'
    averagetype: 'median','mean','mode'
    """
    # PARAMETER SET
    specific = 'MAY18_JUNE18'
    component = '*'
    averagetype = 'median'
    colormap = 'gist_rainbow'
    # //PARAMETER SET


    # various markers and linestyles for easy visualization
    styles=['-']
    markers = ['x','o','D','v','s','*','+']

    # gather data
    stations,periods,averages = multi_component_ppsd(comp=component,
                                                     specific=specific,
                                                     avg=averagetype)


    f,ax = prime_plot()
    color_dict = build_color_dictionary(map='gist_rainbow',
                                        num_of_colors=len(stations))
    for i,(avgs,stas) in enumerate(zip(averages,stations)):
        sta_num = int(stas[2:])
        plt.plot(periods,avgs,
                 c=color_dict[i],
                 markersize=2,
                 markeredgecolor='k',
                 markeredgewidth='.5',
                 marker=markers[i%len(markers)],
                 linestyle=styles[i%len(styles)],
                 label="{} {}".format(stas,station_dict(stas))
                 )

    # final plotting adjustments
    plt.title("{s} HH{c} {a}".format(s=specific,c=component,a=averagetype))
    plt.legend(prop={'size': 4},ncol=4)

    if save:
        figname = "placeholder"
        plt.savefig(figname)
    if show:
        plt.show()

if __name__ == "__main__":
    plot_ppsd()
