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

def setup_plot(number_of_files,twax=True):
    """dynamically set up plots according to number of files
    """
    f = plt.figure(figsize=(11.69,8.27),dpi=75)
    nrows,ncols=number_of_files,1
    height_ratios = [1] * (number_of_files-1)
    height_ratios += [3]
    gs = gridspec.GridSpec(nrows,ncols,height_ratios=height_ratios,hspace=0)

    # create gridspec subplots, sharex with the first axis
    axes,twaxes = [],[]
    for i in range(number_of_files):
        if i == 0:
            ax = plt.subplot(gs[i])
        else:
            ax = plt.subplot(gs[i],sharex=axes[0])
        twax = None
        if twax:
            twinax = ax.twinx()
            twaxes.append(twinax)
        pretty_grids(ax)
        axes.append(ax)

    # remove x tick labels except for last axis
    for ax in axes[0:-1]:
        plt.setp(ax.get_xticklabels(),visible=False)

    return f,axes,twaxes

def create_x_axes(st,Rm=None):
    """create a common x axis for plotting
    """
    start = st[0].stats.starttime.hour
    end = st[0].stats.endtime.hour
    if start > end:
        end += 12
    x_st = np.linspace(start,end,len(st[0].data))
    if Rm:
        x_Rm = np.linspace(start,end,len(Rm))
    else:
        x_Rm = None

    return x_st,x_Rm

def points_over_sigma(Rm):
    """determine which points lie over 2sigma for minute arrays
    """
    meanstdfilepath = pathnames()['data'] + 'TEROR/MEANSTD_night_2338.npz'
    meanstdfile = np.load(meanstdfilepath)
    mean = meanstdfile['Rm_mean']
    std = meanstdfile['Rm_sigma']
    two_sigma = mean + (2*std)

    for i,rmpoint in enumerate(Rm):
        if rmpoint < two_sigma:
            Rm[i] = np.nan

    return Rm

def magnifytremors(date):
    """dynamically plot waveforms and envelopes depending on number of files
    """
    pickle_files,npz_files = collect_files(date)
    NoF = len(pickle_files)
    f,axes,twaxes = setup_plot(NoF+1)
    color_dictionary = build_color_dictionary(num_of_colors=NoF)

    for i,(pick,npz,ax,twax) in enumerate(
                        zip(pickle_files,npz_files,axes[:-1],twaxes[:-1])):
        Rm = parse_npz_file(npz,choice='Rm')
        st,st_envelope = process_data(pick,bounds=[2,8])

        x_st,x_Rm = create_x_axes(st,Rm)
        Rm = points_over_sigma(Rm)
        ax.plot(x_st,st[0].data,c=color_dictionary[i])

        twax.scatter(x_Rm,Rm,marker='o',c='k',s=1.5,)
        station = st[0].get_id().split('.')[1]
        ax.set_ylabel(station)
        if station == "RD07":
            pick_surface = pick
        ax.set_ylim([-np.std(st[0].data)*10,np.std(st[0].data)*10])


    axes[-1].set_xlabel('Time [s]')

    # plot surface wave band of last station plotted
    st_surface,st_envelope = process_data(pick_surface,bounds=[1/30,1/6])
    axes[-1].plot(x_st,st_surface[0].data,c='k')

    plt.show()

def windowtremor(date,window=60*5):
    """plot waveforms in small time windows for more in depth viewing
    window should be given in seconds
    """
    pickle_files,npz_files = collect_files(date)
    NoF = len(pickle_files)
    color_dictionary = build_color_dictionary(num_of_colors=NoF)

    # read in all the data first
    DAT,STA = [],[]
    for i,pick in enumerate(pickle_files):
        st,st_envelope = process_data(pick,bounds=[2,8])
        station = st[0].get_id().split('.')[1]
        DAT.append(st[0].data)
        STA.append(station)
        if station == "RD07":
            st_surface,st_envelope = process_data(pick,bounds=[1/30,1/6])
    window *= st[0].stats.sampling_rate
    x_st,_ = create_x_axes(st)

    # sequential plots for each time window
    for window_start in range(0,len(x_st),int(window)):
        f,axes,_ = setup_plot(NoF+1,twax=False)
        window_end = window_start + window
        ws,we = int(window_start),int(window_end)
        for i,(tr,sta,ax) in enumerate(zip(DAT,STA,axes[:-1])):
            ax.plot(x_st[ws:we],tr[ws:we],c=color_dictionary[i])
            ax.set_ylabel(sta)
        axes[-1].plot(x_st[ws:we],st_surface[0].data[ws:we],c='k')
        axes[-1].set_xlabel('Time [h]')
        plt.show()
        plt.close()

if __name__ == "__main__":
    windowtremor(date='2017-251')
