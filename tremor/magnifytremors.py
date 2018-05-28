"""indepth investigation of tremor detections by pyfreqscan
"""
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from obspy import read, UTCDateTime
from pyfreqscan import waveform_envelope

sys.path.append("../modules")
from getdata import pathnames
from plotmod import pretty_grids, build_color_dictionary
    
def collect_files(date):
    """grab datafile by dates
    """
    date = UTCDateTime(date)
    pickle = pathnames()['data'] + "TEROR/{y}/XX/*/pickle".format(y=date.year)
    pickle_files = glob.glob(os.path.join(pickle,'*{}*'.format(date.julday)))
        
    npz = pathnames()['data'] + "TEROR/{y}/XX/*/npz".format(y=date.year)
    npz_files = glob.glob(os.path.join(npz,'*{}*'.format(date.julday)))
    
    if not pickle_files:
        print("No files found")
        return None
        
    return pickle_files, npz_files
    
def process_data(fid,bounds):
    """read in data and do minor preprocessing. pickle files already have 
    response removed, are preprocessed, and are trimmed to NZ local nighttime 
    if _night tag, so only filtering required
    """    
    st = read(fid)
    st.filter('bandpass',freqmin=bounds[0],freqmax=bounds[1])
    st_envelope = waveform_envelope(st[0].data)
    
    return st, st_envelope
    
def parse_npz_file(fid,choice='Rm'):
    """load in npz file and parse out by request
    """
    npz_file = np.load(fid)
    return npz_file[choice]
        
def setup_plot(number_of_files):
    """dynamically set up plots according to number of files
    """
    f = plt.figure(figsize=(11.69,8.27),dpi=75)
    nrows,ncols=number_of_files,1
    height_ratios = [2] * number_of_files
    gs = gridspec.GridSpec(nrows,ncols,height_ratios=height_ratios,hspace=0)
    
    # create gridspec subplots, sharex with the first axis
    axes,twaxes = [],[]
    for i in range(number_of_files):
        if i == 0:
            ax = plt.subplot(gs[i])
        else:
            ax = plt.subplot(gs[i],sharex=axes[0])
        twinax = ax.twinx()
        pretty_grids(ax)
        axes.append(ax)
        twaxes.append(twinax)

    # remove x tick labels except for last axis
    for ax in axes[0:-1]:
        plt.setp(ax.get_xticklabels(),visible=False) 
    
    return f,axes,twaxes

def create_x_axes(st,Rm):
    """create a common x axis for plotting
    """
    start = st[0].stats.starttime.hour
    end = st[0].stats.endtime.hour
    if start > end:
        end += 12
    x_st = np.linspace(start,end,len(st[0].data))
    x_Rm = np.linspace(start,end,len(Rm))
    
    return x_st,x_Rm
    
def magnifytremors(date):
    """dynamically plot waveforms and envelopes depending on number of files
    """
    pickle_files,npz_files = collect_files(date)
    NoF = len(pickle_files)
    f,axes,twaxes = setup_plot(NoF)
    color_dictionary = build_color_dictionary(num_of_colors=NoF)
    
    for i,(pick,npz,ax,twax) in enumerate(
                                    zip(pickle_files,npz_files,axes,twaxes)):
        Rm = parse_npz_file(npz,choice='Rm')
        st,st_envelope = process_data(pick,bounds=[2,8])
        x_st,x_Rm = create_x_axes(st,Rm)
        ax.plot(x_st,st[0].data,c=color_dictionary[i])
        twax.plot(x_Rm,Rm,c='k')
        
    plt.show()
        

        
if __name__ == "__main__":
    magnifytremors(date='2017-251')
        
    
                                                                  
    