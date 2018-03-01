"""24/1/18
Use geonet permanent station data for noise spectra analysis
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
from getdata import geonet_internal, fdsn_download, pathnames

# user input arguments
parser = argparse.ArgumentParser(description='Download create PPSD for \
 station and timeframe $ python noise_analysis.py --station TBAS --channel\
 HHZ --start 2015-01-01 --end 2016-01-01')
parser.add_argument('--station', help='Instrument station, default = KNZ',
                    type=str,default='KNZ')
parser.add_argument('--channel', help='Instrument channel, default = HHZ,\
                    other options, BNE, EHZ',
                    type=str, default='Z')
parser.add_argument('--start', help='Starttime YYYY-MM-DD default = \
                    2015-01-01',type=str, default='2015-01-01')
parser.add_argument('--end', help='Endtime, default = 2015-01-01',type=str,
                    default='2015-01-02')
parser.add_argument('--dec', help='Decimate trace, default 0',type=int,
                    default=0)

# parse arguments
arg = parser.parse_args()
station = arg.station.upper()
channel = arg.channel.upper()
start = UTCDateTime(arg.start)
end = UTCDateTime(arg.end)
decimateby = arg.dec

print("\nTimeframe set: {} through {}".format(start.date,end.date))

# set instrument choice on command line
network = 'NZ'
if station[:2].upper() == 'RD':
    network = 'XX'

instrument_id = '{net}.{sta}.*.{cha}'.format(net=network,
                                            sta=station,
                                            cha=channel)
print("Analyzing noise spectra for {}".format(instrument_id))

# filepaths of geonet archive
if network == 'NZ':
    if pathnames()["where"] == "GNS":
        data_files, resp_file = geonet_internal(station=station,
                                            channel=channel,
                                            start=start,
                                            end=end)

# RDF temporary network - includes all files, does not filter by time
elif network == 'XX':
    d1 = pathnames()["rdf"] + "July2017_Sep2017/DATA_ALL/"
    d2 = pathnames()["rdf"] + "Sep2017_Nov2017/DATA_ALL"
    d3 = pathnames()["rdf"] + "Nov2017_Jan2018/DATA_ALL/"
    data_files = []
    for D in [d1,d2,d3]:
        data_files += glob.glob(D + '{sta}*{cha}*'.format(sta=station,
                                                        cha=channel))
    data_files.sort()
    print("++ {} data files found".format(len(data_files)))
    resp_file = pathnames()["rdf"] + "RDF_Array/DATALESS.RDF.XX"

# read in response file, set decimate parameter
print("Response: {}".format(resp_file[0]))
inv = read_inventory(resp_file[0])

# initialize PPSD with first datafile
print("1/{} Initializing with data file: ".format(len(data_files)),
                                                    data_files[0],end='... ')
start = time.time()
st = read(data_files[0])
if decimateby != 0:
    st.decimate(decimateby)
tr = st.select(channel=channel)[0]
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

output_filename = '{sta}.{cha}.{yearS}.{day_start:0>3}-{yearE}.{day_end:0>3}'.format(
                                                        sta=station,
                                                        cha=channel,
                                                        yearS=year_start,
                                                        yearE=year_end,
                                                        day_start=jday_start,
                                                        day_end=jday_end)

npz_filepath = './ppsd_arrays/{}.npz'.format(output_filename)
png_filepath = './output_plots/ppsd_plots/{}.png'.format(output_filename)

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

# ignore overwrite protection
ppsd.save_npz(npz_filepath)
print("saved .npz file: ",npz_filepath)

# ppsd.plot(filename=png_filepath,cmap=pqlx,show_mean=True)
# print("saved .png file: ",png_filepath)
