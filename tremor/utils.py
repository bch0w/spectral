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

    !!! base code copied and modified from obspy.imaging.waveform.__plot_min_max
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
    pixel_length = pixel_length
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
    net,sta,loc,cha,year,jday = code_set.split('.')
    outpath = pathnames()['data'] + 'TEROR'
    nightcheck = ""
    if night:
        nightcheck = "_night"
    outfile = "{n}.{s}.{l}.{y}.{j}.TEORRm{N}.{f}".format(n=net,
                                                      s=sta,
                                                      l=loc,
                                                      y=year,
                                                      j=jday,
                                                      N=nightcheck,
                                                      f='{f}')
    output = os.path.join(outpath,outfile)
    pickle_path = output.format(f='pickle')
    npz_path = output.format(f='npz')

    # save arrays
    if TEORRm:
        T,E,O,R,Rm,Rh = TEORRm
        st.write(pickle_path,format="PICKLE")
        np.savez(npz_path,T=T,E=E,O=O,R=R,Rm=Rm,Rh=Rh)
        return True

    # check arrays
    else:
        if (os.path.exists(npz_path) and os.path.exists(pickle_path)):
            return {"npz":npz_path,"pickle":pickle_path}
        else:
            return False

def already_processed():
    """check what days have already been processed
    """
    import glob
    output_path = pathnames()['data'] + 'TEROR/*npz'
    allfiles = glob.glob(output_path)
    jday_list,sta_list = [],[]
    for fid in allfiles:
        fid = os.path.basename(fid)
        net,sta,loc,year,jday,_,_ = fid.split('.')
        jday_list.append(jday)
        sta_list.append(sta)

    jday_list,sta_list = zip(*sorted(zip(jday_list,sta_list), key=lambda x:x[0]))
    for j,s in zip(jday_list,sta_list):
        print(j,s)
    sys.exit()
