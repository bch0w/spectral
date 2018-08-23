"""for statistical plotting of data, such as histograms and bar charts
"""
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import sys
sys.path.append('..')
import corkBoard
sys.path.append('../../modules')
from procmod import myround


class Depicter:
    """a class used to plot information from Tacks
    """
    def __init__(self,**kwargs):
        self.tack = kwargs.get('tack')

    def _set_params(self):
        import matplotlib as mpl
        mpl.rcParams['font.size'] = 12
        mpl.rcParams['lines.linewidth'] = 1.25
        mpl.rcParams['lines.markersize'] = 10
        mpl.rcParams['axes.linewidth'] = 2

    def plot_misfit_histogram(self,*args,**kwargs):
        """create a histogram from the misfit values
        """
        # only import plotting functions if necessary
        self._set_params()
        
        binsize = kwargs.get('binsize',0.1)
        m0 = kwargs.get('m0',None)
        m_a = kwargs.get('m_a',None)

        # gather misfit values
        misfits = np.fromiter(self.tack.misfit_values.values(),dtype="float")
        maxmisfit = myround(misfits.max(),base=1,choice="up")
        n,bins,patches = plt.hist(x=misfits,
                                  bins=len(np.arange(0,maxmisfit,binsize)),
                                  range=(0,maxmisfit),
                                  color="orange",
                                  histtype="bar",
                                  edgecolor="black",
                                  linewidth=1.5,
                                  zorder=10)

        # plot attributes
        plt.xlabel("Misfit value")
        plt.ylabel("Count (N={})".format(len(misfits)))
        plt.title("{eid}$_{mod}$ Misfits ".format(eid=self.tack.id,mod="{m00}"))
        plt.grid(linewidth=1.0,which='both',zorder=1)
        plt.xlim([-0.05,bins.max()+0.05])
        plt.ylim([0,max(n)+0.5])
        plt.show()
        
    def plot_cc_time_shift_histogram(self,*args,**kwargs):
        """create a histogram of cross correlation time shifts
        """
        self._set_params()
        binsize = kwargs.get('binsize',0.1)
        
        # gather cc time shift values
        cc_time_shift = np.fromiter(self.tack.cc_time_shifts.values(),
                                                                dtype="float")
        max_cc_time_shift = myround(
                    max([abs(cc_time_shift.max()),abs(cc_time_shift.min())]),
                    base=0.1,choice="up")
        
        n,bins,patches = plt.hist(x=cc_time_shift,
                                  bins=len(
                                        np.arange(0,max_cc_time_shift,binsize)),
                                  range=(-max_cc_time_shift,max_cc_time_shift),
                                  color="orange",
                                  histtype="bar",
                                  edgecolor="black",
                                  linewidth=1.5,
                                  zorder=10)

        # annotate mean and one sigma, plot it
        mean = np.mean(cc_time_shift)
        onesigma = np.std(cc_time_shift)
        mean_one_sigma = "$\mu$ + $\sigma$ = {m:.2f} $\pm$ {s:.2f}s".format(
                                                    m=mean,
                                                    s=onesigma)
        plt.text(x=-max_cc_time_shift+.5,y=n.max()-1,s=mean_one_sigma,
                 bbox=dict(facecolor='w',alpha=0.5),zorder=11)
        plt.axvline(x=mean,ymin=n.min(),ymax=n.max(),linestyle='-',color='k')
        for i in [-1,1]:
            plt.axvline(x=i*onesigma,ymin=n.min(),ymax=n.max(),linestyle='--',
                        color='k')
        
        # plot attributes
        plt.xlabel("CC Time Shift [s]")
        plt.ylabel("Count (N={})".format(len(cc_time_shift)))
        plt.title("{eid}$_{mod}$ CC Time Shifts ".format(eid=self.tack.id,
                                                           mod="{m00}"))
        plt.grid(linewidth=1.0,which='both',zorder=1)
        plt.xlim([-(max_cc_time_shift+0.05),max_cc_time_shift+0.05])
        plt.ylim([0,max(n)+0.5])
        plt.show()
