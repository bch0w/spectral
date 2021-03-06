"""utility function for pyfreqscan
"""
import os
import sys
import numpy as np
from matplotlib.dates import date2num

sys.path.append("../modules")
from getdata import pathnames

def create_min_max(tr,pixel_length=10):
    """
    Creates new data using a min/max approach that calculated the minimum and
    maximum values of each "pixel". Much faster plotting of large datasets.

    !!! base code copied and modified 
    !!! from obspy.imaging.waveform.__plot_min_max
    """
    # Some variables to help calculate the values.
    starttime = date2num(tr.stats.starttime.datetime)
    endtime = date2num(tr.stats.endtime.datetime)
    # The same trace will always have the same sampling_rate.
    sampling_rate = tr.stats.sampling_rate
    # width of x axis in seconds
    x_width = endtime - starttime
    # number of samples that get represented by one min-max pair
    # width = 800 # guessing
    # pixel_length = int(
    #     np.ceil((x_width * sampling_rate + 1) / width))
    # Loop over all the traces. Do not merge them as there are many samples
    # and therefore merging would be slow.
    trace_length = len(tr.data)
    pixel_count = int(trace_length // pixel_length)
    remaining_samples = int(trace_length % pixel_length)
    remaining_seconds = remaining_samples / sampling_rate
    # Reference to new data array which does not copy data but can be
    # reshaped.
    if remaining_samples:
        data = tr.data[:-remaining_samples]
    else:
        data = tr.data
    data = data.reshape(pixel_count, pixel_length)
    min_ = data.min(axis=1) * tr.stats.calib
    max_ = data.max(axis=1) * tr.stats.calib
    # Calculate extreme_values and put them into new array.
    if remaining_samples:
        extreme_values = np.empty((pixel_count + 1, 2), dtype=np.float)
        extreme_values[:-1, 0] = min_
        extreme_values[:-1, 1] = max_
        extreme_values[-1, 0] = \
            tr.data[-remaining_samples:].min() * tr.stats.calib
        extreme_values[-1, 1] = \
            tr.data[-remaining_samples:].max() * tr.stats.calib
    else:
        extreme_values = np.empty((pixel_count, 2), dtype=np.float)
        extreme_values[:, 0] = min_
        extreme_values[:, 1] = max_
    # Finally plot the data.
    start = date2num(tr.stats.starttime.datetime)
    end = date2num(tr.stats.endtime.datetime)
    if remaining_samples:
        # the last minmax pair is inconsistent regarding x-spacing
        x_values = np.linspace(start, end - remaining_seconds,
                               num=extreme_values.shape[0] - 1)
        x_values = np.concatenate([x_values, [end]])
    else:
        x_values = np.linspace(start, end, num=extreme_values.shape[0])
    x_values = np.repeat(x_values, 2)
    y_values = extreme_values.flatten()

    return x_values, y_values

def __z2nan(array):
    """convert zeros in an array to nan for plotting use
    """
    array[array==0] = np.nan

    return array

def check_save(code_set,st=None,TEORRm=None,night=False):
    """either check if processing has been run before, or save files to
    specific path. [Check] initiated by default, [save] initiated if the
    function is given an argument for TEORRm
    :type code_set: str
    :param code_set: instrument code, set in main
    :type st: obspy stream
    :param st: stream containing preprocessed data
    :type TEORRm: list of numpy arrays
    :param TEORRm: list of arrays containing filtered waveform data
    """
    # set up pathing
    net,sta,loc,cha,d,year,jday = code_set.split('.')
    outpath = os.path.join(pathnames()['data'],'TEROR',year,net,sta,'{f}')
    nightcheck = ""
    if night:
        nightcheck = "_night"
    outfile = "{n}.{s}.{l}.{y}.{j}{N}.{f}".format(n=net,
                                                      s=sta,
                                                      l=loc,
                                                      y=year,
                                                      j=jday,
                                                      N=nightcheck,
                                                      f='{f}')
    output = os.path.join(outpath,outfile)
    pickle_path = output.format(f='pickle')
    npz_path = output.format(f='npz')
    for pathcheck in [pickle_path,npz_path]:
        if not os.path.exists(os.path.dirname(pathcheck)):
            os.makedirs(os.path.dirname(pathcheck))

    # save stream as a pickle to preserve non-mseed components
    # save passband arrays into npz file, dynamically sort out dict components**
    if st:
        if not os.path.exists(pickle_path):
            st.write(pickle_path,format="PICKLE")
        if TEORRm:
            np.savez(npz_path,**TEORRm)
            return True
        else:
            return False

    # check arrays
    else:
        if not os.path.exists(npz_path):
            npz_path = False
        if not os.path.exists(pickle_path):
            pickle_path = False
        return {"npz":npz_path,"pickle":pickle_path}

def already_processed():
    """check what days have already been processed, pretty print by day number
    """
    import glob
    output_path = os.path.join(pathnames()['data'],
                                    'TEROR','201?','*','*','pickle','*.pickle')
    allfiles = glob.glob(output_path)
    jday_list,sta_list = np.array([]),np.array([])
    for fid in allfiles:
        fid = os.path.basename(fid)
        net,sta,loc,year,jday,_ = fid.split('.')
        jday_list = np.append(jday_list,jday)
        sta_list = np.append(sta_list,sta)

    jday_set = sorted(set(jday_list))
    for J in jday_set:
        indices = np.where(jday_list==J)[0]
        print(J,end=": ")
        for i in indices:
            print(sta_list[i],end=" ")
        print(' ')
    sys.exit()
