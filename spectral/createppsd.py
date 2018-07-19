"""24/1/18
Use geonet permanent sta data for noise spectra analysis
in accordance with McNamara2004 using the obspy PPSD package
"""
import os
import sys
import glob
import time
import argparse
from obspy import read
from obspy import read_inventory
from obspy.signal import PPSD
from obspy import UTCDateTime
from obspy.imaging.cm import pqlx
import matplotlib
matplotlib.use('agg')

# internal modules
sys.path.append('../modules/')
import getdata
from getdata import pathnames


def gather_data(sta,cha,start,end):
    """collect data and response information, currently only getting RDF data
    """
    mseed_files, response_file = getdata.get_fathom(station=sta,channel=cha,
                                                        start=start,end=end)

    return mseed_files, response_file

def analyze_noise(data_files,response,decimateby=5):
    """run through data files and create PPSD objects using obsy
    """
    data_files.sort()
    print("++ {} data files".format(len(data_files)))
    inv = read_inventory(response)
    # initialize PPSD with first datafile
    print("1/{} Initializing with data file: ".format(len(data_files)),
                                    os.path.basename(data_files[0]),end='... ')
    start = time.time()
    st = read(data_files[0])
    if decimateby != 0:
        st.decimate(decimateby)
    ppsd = PPSD(st[0].stats, metadata=inv)
    ppsd.add(st)
    year_start = st[0].stats.starttime.year
    jday_start = st[0].stats.starttime.julday
    end = time.time()
    print("complete ({}s)".format(round(end-start,2)))

    # loop over rest of datafiles and add to ppsd
    for i,filename in enumerate(data_files[1:]):
        print('{0}/{1} {2}'.format(i+2,len(data_files),
                                        os.path.basename(filename)),end='... ')
        try:
            start = time.time()
            st = read(filename)
            if decimateby != 0:
                st.decimate(decimateby)
            ppsd.add(st)
            end = time.time()
            print("complete ({}s)".format(round(end-start,2)))
        except Exception as e:
            print(e)
            pass

    return ppsd

def single_process():
    start = '2018-028'
    end = '2018-057'
    station = 'RD01'
    channel = 'HHZ'
    decimate_by = 5
    data,response = gather_data(sta=station,cha=channel,start=start,end=end)
    ppsd = analyze_noise(data,response,decimate_by)

    # save ppsd and figure
    output_name = '{s}.{c}.{st}.{en}.db{d}'.format(s=station,
                                                   c=channel,
                                                   st=start,
                                                   en=end,
                                                   d=decimate_by)

    npz_filepath = pathnames()['ppsd'] + '{}.npz'.format(output_name)
    png_filepath = pathnames()['plots'] + \
                    'spectral/ppsd_plots/PPSD_HOLD/{}.png'.format(output_name)


    ppsd.plot(filename=png_filepath,cmap=pqlx,show_mean=True,show=True)

def mass_process():
    """main processing function
    """
    # common parameter set
    start = '2017-199'
    end = '2018-001'
    decimate_by = 5
    for i in range(1,23):
        # if i in [10,11,17]:
        #     continue
        station = 'RD{:0>2}'.format(i)
        # for channel in ['HHZ','HHN','HHE']:
        for channel in ['HHN','HHE']:
            data,response = gather_data(sta=station,cha=channel,
                                                        start=start,end=end)
            if not data:
                continue
            ppsd = analyze_noise(data,response)

            # save ppsd and figure
            output_name = '{s}.{c}.{st}.{en}.db{d}'.format(s=station,
                                                           c=channel,
                                                           st=start,
                                                           en=end,
                                                           d=decimate_by)

            npz_filepath = pathnames()['ppsd'] + '{}.npz'.format(output_name)
            png_filepath = pathnames()['plots'] + \
                            'spectral/ppsd_plots/PPSD_HOLD/{}.png'.format(
                                                                    output_name)
            ppsd.save_npz(npz_filepath)

if __name__ == '__main__':
    mass_process()
