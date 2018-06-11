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

def station_dict(code):
    station_dict = {"RD01":"PRWZ", "RD02":"ANWZ", "RD03":"TURI", "RD04":"PORA",
                    "RD05":"MNHR", "RD06":"DNVZ", "RD07":"WPAW", "RD08":"RAKW",
                    "RD09":"MCHZ", "RD10":"CKHZ", "RD11":"KAHU", "RD12":"KWHZ",
                    "RD13":"KERE", "RD14":"PNUI", "RD15":"WPUK", "RD16":"OROA",
                    "RD17":"TEAC", "RD18":"RANC","RD19":"MATT","RD20":"KAHU2",
                    "RD21":"TEAC2"}
    return station_dict[code]

def gather_ppsd(specific,comp='Z',avg='median'):
    """grab folders of data. comp can be wildcard or N,E or Z
    """
    npz_path = pathnames()['ppsd']
    filename = '*HH{}*.npz'.format(comp)
    files = glob.glob(os.path.join(npz_path,specific,filename))
    files.sort()

    averages,stations = [],[]
    for fid in files:
        sta,cha,date1,date2,db,fmt = os.path.basename(fid).split(".")
        ppsd = PPSD.load_npz(fid)
        if avg == 'mean':
            average = ppsd.get_mean()
        elif avg == 'mode':
            average = ppsd.get_mode()
        elif avg  == 'median':
            average = ppsd.get_percentile(percentile=50)
        averages.append(average[1])
        stations.append(sta)

    periods = average[0]

    return stations,periods,averages

def compare_ppsd(comp='Z',avg='median'):
    """grab folders of data. comp can be wildcard or N,E or Z
    """
    npz_path = pathnames()['ppsd']
    specific1 = 'RDF_JanMar_deployment/'
    specific2 = 'RDF_MarMay_deployment/'
    specific_list = [specific1,specific2]
    averages,stations = [],[]
    for i,specific in enumerate(specific_list):
        filename = '*HH{}*.npz'.format(comp)
        files = glob.glob(os.path.join(npz_path,specific,filename))
        files.sort()

        for fid in files:
            sta,cha,date1,date2,db,fmt = os.path.basename(fid).split(".")
            ppsd = PPSD.load_npz(fid)
            if avg == 'mean':
                average = ppsd.get_mean()
            elif avg == 'mode':
                average = ppsd.get_mode()
            elif avg  == 'median':
                average = ppsd.get_percentile(percentile=50)
            averages.append(average[1])
            stations.append(sta)
        periods = average[0]

    return stations,periods,averages

def prime_plot():
    """set up plot objects
    """
    f = plt.figure(figsize=(9,5),dpi=200)
    ax = f.add_subplot(111)
    pretty_grids(ax)

    # plot lines for noise models and microseisms
    nlnm_x,nlnm_y = get_nlnm()
    nhnm_x,nhnm_y = get_nhnm()
    plt.plot(nhnm_x,nhnm_y,'gray',alpha=0.4)
    plt.plot(nlnm_x,nlnm_y,'gray',alpha=0.4)
    ax.fill_between(nhnm_x,nlnm_y,nhnm_y,facecolor='gray',alpha=0.2)

    plt.xlim([0.2,100])
    plt.ylim([nlnm_y.min(),-90])
    plt.xscale("log")
    plt.xlabel("Period [s]")
    plt.ylabel("Amplitude [m^2/s^4/Hz][dB]")
    plt.grid()

    return f, ax


def plot_ppsd(save=False,show=True):
    """plotting function
    """
    styles = ['-.','--','-']
    color_dict = build_color_dictionary(map='viridis',num_of_colors=21)

    # compare stations
    stations,periods,averages = compare_ppsd()
    f,ax = prime_plot()

    for i,(avgs,stas) in enumerate(zip(averages,stations)):
        sta_num = int(stas[2:])
        plt.plot(periods,avgs,
                 label="{} {}".format(stas,station_dict(stas)))

    if save:
        figname = "placeholder"
        plt.savefig(figname)
    if show:
        plt.show()

def plot_by_station():
    color_dict = build_color_dictionary(map='viridis',num_of_colors=21)

    # compare stations
    specific1 = 'RDF_JanMar_deployment'
    specific2 = 'RDF_MarMay_deployment'
    S1,periods,A1 = gather_ppsd(specific1)
    S2,periods,A2 = gather_ppsd(specific2)

    for i in range(1,22):
        station_number = 'RD{:0>2}'.format(i)
        f,ax = prime_plot()
        for S,A,specific,L in zip([S1,S2],[A1,A2],[specific1,specific2],['-','--']):
            try:
                index = S.index(station_number)
                plt.plot(periods,A[index],label=specific,linestyle=L)
                plt.legend()
                plt.title(station_number)
            except ValueError:
                continue
        plt.show()



if __name__ == "__main__":
    plot_by_station()
