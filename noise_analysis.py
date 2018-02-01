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
from getdata_fdsn import geonet_data

# user input arguments
parser = argparse.ArgumentParser(description='Download create PPSD for \
 station and timeframe $ python noise_analysis.py --station TBAS -comp\
 z --start 2015-01-01 --end 2016-01-01')
parser.add_argument('--station', help='Instrument station, default = KNZ',
                    type=str,default='KNZ')
parser.add_argument('--comp', help='Instrument component, default = Z',
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
comp = arg.comp.upper()
start = UTCDateTime(arg.start)
end = UTCDateTime(arg.end)
decimateby = arg.dec

print("\nTimeframe set: {} through {}".format(start.date,end.date))

# set instrument choice on command line
if station[-1].upper() == 'Z':
    channel = 'HH' + comp
elif station[-1].upper() == 'S':
    channel = 'BN' + comp
instrument_id = 'NZ.{}.*.{}'.format(station,channel)
print("Analyzing noise spectra for {}".format(instrument_id))

# filepaths of geonet archive
data_files, resp_file = geonet_data(station=station,
                                    comp=channel[-1],
                                    start=start,
                                    end=end)

# read in response file, set decimate parameter
inv = read_inventory(resp_file)

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
year_start,jday_start = os.path.basename(data_files[0]).split('.')[-2:]
year_end,jday_end = os.path.basename(data_files[-1]).split('.')[-2:]

output_filename = '{sta}.{cha}.{year}.{day_start}-{day_end}'.format(sta=station,
                                                        cha=channel,
                                                        year=year_start,
                                                        day_start=jday_start,
                                                        day_end=jday_end)

npz_filepath = './ppsd_arrays/{}.npz'.format(output_filename)
png_filepath = './figures/ppsd_plots/{}.png'.format(output_filename)

# check npz file exists
if os.path.exists(npz_filepath):
    npz_check = input('npz path exists, overwrite? [y]/n')
    if npz_check == 'n':
        ppsd.save_npz(npz_filepath)
else:
    ppsd.save_npz(npz_filepath)
    print("saved .npz file: ",npz_filepath)

# check png file exists
if os.path.exists(png_filepath):
    png_check = input('png path exists, overwrite? [y]/n')
    if png_check == 'n':
        import pdb;pdb.set_trace()
        sys.exit('Aborted')

ppsd.plot(filename=png_filepath,cmap=pqlx,show_mean=True)
print("saved .png file: ",png_filepath)
