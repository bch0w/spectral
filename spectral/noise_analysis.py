"""24/1/18
Use geonet permanent sta data for noise spectra analysis
in accordance with McNamara2004 using the obspy PPSD package
"""
import os
import sys
sys.path.append('../modules/')
import glob
import time
import argparse
from obspy import read
from obspy import read_inventory
from obspy.signal import PPSD
from obspy import UTCDateTime
from obspy.imaging.cm import pqlx

import getdata
from getdata import pathnames


# user input arguments
parser = argparse.ArgumentParser(description='Download create PPSD for \
 sta and timeframe $ python noise_analysis.py --sta TBAS --cha\
 HHZ --start 2015-01-01 --end 2016-01-01')
parser.add_argument('--code', help='Instrument Code, i.e. XX.RD01..HHZ or\
                        NZ.PUZ..HHE', type=str, default='NZ.PUZ..HHE')
parser.add_argument('--start', help='Starttime YYYY-MM-DD default = \
                    2015-01-01',type=str, default='2015-01-01')
parser.add_argument('--end', help='Endtime, default = 2015-01-01',type=str,
                    default='2015-01-02')
parser.add_argument('--dec', help='Decimate trace, default 0',type=int,
                    default=0)

# parse arguments
arg = parser.parse_args()
code = arg.code
start = UTCDateTime(arg.start)
end = UTCDateTime(arg.end)
decimateby = arg.dec

net,sta,loc,cha = code.split('.')
print("Analyzing noise spectra for {}".format(code))

# filepaths of geonet archive
if net == 'NZ':
    if pathnames()["where"] == "GNS":
        print("\nTimeframe set: {} through {}".format(start.date,end.date))
        data_files, resp_file = geonet_internal(station=sta,
                                            channel=cha,
                                            start=start,
                                            end=end)

# RDF temporary net - includes all files, does not filter by time
elif net == 'XX':
    d1 = pathnames()['rdf'] + "July2017_Sep2017/DATA_ALL/",
    d2 = pathnames()['rdf'] + "Sep2017_Nov2017/DATA_ALL/",
    d3 = pathnames()['rdf'] + "Nov2017_Jan2018/DATA_ALL/",
    d4 = pathnames()['rdf'] + "Jan2018_Mar2018/DATA_ALL/"
    array_list = [d4]
    data_files = []
    for D in array_list:
        data_files += glob.glob(D + '{sta}*{cha}*'.format(sta=sta,
                                                        cha=cha))
    data_files.sort()
    print("++ {} data files found".format(len(data_files)))
    resp_file = pathnames()['rdf'] + "DATALESS.RDF.XX"

# read in response file, set decimate parameter
print("Response: {}".format(resp_file))
inv = read_inventory(resp_file)

# initialize PPSD with first datafile
print("1/{} Initializing with data file: ".format(len(data_files)),
                                                    data_files[0],end='... ')
start = time.time()
import ipdb;ipdb.set_trace()
st = read(data_files[0])
if decimateby != 0:
    st.decimate(decimateby)
tr = st.select(channel=cha)[0]
ppsd = PPSD(tr.stats, metadata=inv)
ppsd.add(st)
year_start = st[0].stats.starttime.year
jday_start = st[0].stats.starttime.julday
end = time.time()
print("complete ({}s)".format(round(end-start,2)))

# loop over rest of datafiles and add to ppsd
for i,filename in enumerate(data_files[1:]):
    print('{0}/{1} {2}'.format(i+2,len(data_files),filename),end='... ')
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

# set output name for saving .npz and f# save data as .npz array, save plot as .png
# i'm assuming here the year stays the same
year_end = st[0].stats.starttime.year
jday_end = st[0].stats.starttime.julday

output_template = '{sta}.{cha}.{yearS}.{day_start:0>3}-{yearE}.{day_end:0>3}'
output_name = output_template.format(sta=sta,
                                        cha=cha,
                                        yearS=year_start,
                                        yearE=year_end,
                                        day_start=jday_start,
                                        day_end=jday_end)

npz_filepath = pathnames()['ppsd'] + '{}.npz'.format(output_name)
png_filepath = pathnames()['plots'] + 'ppsd_plots/{}.png'.format(output_name)

# ignore overwrite protection
ppsd.save_npz(npz_filepath)
print("saved .npz file: ",npz_filepath)
ppsd.plot(filename=png_filepath,cmap=pqlx,show_mean=True)
print("saved .png file: ",png_filepath)

# overwrite protection
# check npz file exists
# if os.path.exists(npz_filepath):
#     npz_check = input('npz path exists, overwrite? [y]/n')
#     if npz_check == 'n':
#         ppsd.save_npz(npz_filepath)
# else:
#     ppsd.save_npz(npz_filepath)
#     print("saved .npz file: ",npz_filepath)

# check png file exists
# if os.path.exists(png_filepath):
#     png_check = input('png path exists, overwrite? [y]/n')
#     if png_check == 'n':
#         import pdb;pdb.set_trace()
#         sys.exit('Aborted')
