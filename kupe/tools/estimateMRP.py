"""quick plot script for estimating Minimum Resolveable Period for a coarse
(NGLL5) vs fine (NGLL7) mesh
"""
import os
import sys
sys.path.append('../../modules')
from obspy import read, Stream
import numpy as np

from getdata import pathnames
from plotmod import pretty_grids

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1

def estimate_MRP(sta):
    """the function
    """
    basepath = pathnames()['syns'] + '2018p130600_ngllTest'
    objective_filename = "NZ.{s}.{m}XZ.semv.mseed"

    # read in stream objects and mark
    st = Stream()
    for NGLL,marker in zip(['NGLL5','NGLL7'],['B','H']):
        filename = os.path.join(basepath,NGLL,
                                    objective_filename.format(s=sta,m=marker))
        st_tmp = read(filename)
        st_tmp[0].stats.location = NGLL
        st += st_tmp

    # determine common endtime
    commonend = min([st[0].stats.endtime,st[1].stats.endtime])
    commonSR = min([st[0].stats.sampling_rate,st[1].stats.sampling_rate)
    for tr in st:
        tr.trim(tr.stats.starttime,commonend)
        tr.resample(commonSR)
        
    # setup plot
    f,ax = plt.subplots()
    pretty_grids(ax)

    # setup time axis
    t_ngll5 = np.linspace(
                    0,st[0].stats.endtime-st[0].stats.starttime,len(st[0].data))
    t_ngll7 = np.linspace(
                    0,st[1].stats.endtime-st[1].stats.starttime,len(st[1].data))




    step = 0
    st.filter('bandpass',freqmin=1/100,freqmax=1)
    for lowpass in np.arange(0.5,8,0.5):
        st_filter = st.copy()
        st_filter.filter('lowpass',freq=lowpass)
        for tr in st_filter:
            tr.data /= tr.data.max()
            tr.data += step

        plt.plot(t_ngll5,st_filter.select(location='NGLL5')[0].data,'k',label='NGLL5')
        plt.plot(t_ngll7,st_filter.select(location='NGLL7')[0].data,'r',label='NGLL7')
        plt.annotate('lowpass @ {}'.format(lowpass),xy=(0,step))

        step += 1.0

    plt.show()

if __name__ == "__main__":
    estimate_MRP(sta='BFZ')