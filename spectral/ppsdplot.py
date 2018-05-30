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
    
def gather_ppsd(comp='Z',avg='median'):
    """grab folders of data. comp can be wildcard or N,E or Z
    """
    npz_path = pathnames()['ppsd']
    # specific = 'FATHOM/Jan18_Mar18'
    specific = 'FATHOM/Jul17_Jan18'
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
    specific1 = 'Jan18_Mar18'
    specific2 = 'Jul17_Jan18'
    averages,stations = [],[]
    for i,specific in enumerate([specific1,specific2]):
        filename = '*HH{}*.npz'.format(comp)
        files = glob.glob(os.path.join(npz_path,'FATHOM',specific,filename))
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
    
def plot_ppsd(save=False,show=True):
    """plotting function
    """
    f = plt.figure(figsize=(9,5),dpi=200)
    ax = f.add_subplot(111)
    pretty_grids(ax)
    
    styles = ['-.','--','-']
    color_dict = build_color_dictionary(map='viridis',num_of_colors=21)
    
    # plot lines for noise models and microseisms
    nlnm_x,nlnm_y = get_nlnm()
    nhnm_x,nhnm_y = get_nhnm()
    plt.plot(nhnm_x,nhnm_y,'gray',alpha=0.4)
    plt.plot(nlnm_x,nlnm_y,'gray',alpha=0.4)
    ax.fill_between(nhnm_x,nlnm_y,nhnm_y,facecolor='gray',alpha=0.2)
    
    # compare stations
    stations,periods,averages = compare_ppsd()
    stationcheck = [10,11,17,19,20,21]
    stylecheck = ['-','--','-.','-','--','-.']
    colorcheck = ['r','r','r','k','k','k']
    y=0
    for i,(avgs,stas) in enumerate(zip(averages,stations)):
        sta_num = int(stas[2:])
        if sta_num in stationcheck:
            plt.plot(periods,avgs,c=colorcheck[y],
                     linestyle=stylecheck[y],
                     label="{} {}".format(stas,station_dict(stas)))
            y+=1
    
    # plot single timeframe
    # stations,periods,averages = gather_ppsd()
    # for i,(avgs,stas) in enumerate(zip(averages,stations)):
    #     sta_num = int(stas[2:])
    #     if sta_num in stationcheck:
    #         plt.plot(periods,avgs,c=color_dict[sta_num],
    #                  linestyle=styles[i%len(styles)],
    #                  label="{}".format(stas))

    plt.xlim([0.2,100])
    plt.ylim([nlnm_y.min(),-90])
    plt.xscale("log")
    plt.xlabel("Period [s]")
    plt.ylabel("Amplitude [m^2/s^4/Hz][dB]")
    plt.legend(ncol=3,prop={"size":5})
    plt.title("Station Relocation")
    plt.grid()
    if save:
        figname = "placeholder"
        plt.savefig(figname)
    if show:
        plt.show()
    
if __name__ == "__main__":
    plot_ppsd()




