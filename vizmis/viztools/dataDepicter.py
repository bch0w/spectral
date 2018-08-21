"""for statistical plotting of data, such as histograms and bar charts
"""
import numpy as np
import matplotlib.cm as cm

import sys
sys.path.append('..')
import corkBoard
sys.path.append('../../modules')
from procmod import myround


class Depicter:
    """a class used to plot information from Corks
    """
    def __init__(self,**kwargs):
        self.cork = kwargs.get('cork')

    def _set_params(self):
        import matplotlib as mpl
        mpl.rcParams['font.size'] = 12
        mpl.rcParams['lines.linewidth'] = 1.25
        mpl.rcParams['lines.markersize'] = 10
        mpl.rcParams['axes.linewidth'] = 2

    def plot_misfit_histogram(self,*args,**kwargs):
        """create a histogram from the misfit values
        """
        if not self.cork.misfit_values:
            self.cork.collect_misfits()

        # only import plotting functions if necessary
        import matplotlib.pyplot as plt
        self._set_params()
        
        binsize = kwargs.get('binsize',0.1)
        m0 = kwargs.get('m0',None)
        m_a = kwargs.get('m_a',None)

        # gather misfit values
        misfits = np.fromiter(self.cork.misfit_values.values(),dtype="float")
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
        plt.ylabel("Count")
        plt.title("{mod} Misfits ".format(mod="m00"))
        plt.grid(linewidth=1.0,which='both',zorder=1)
        plt.xlim([-0.05,bins.max()+0.05])
        plt.ylim([0,max(n)+0.5])
        plt.show()
