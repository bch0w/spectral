"""24/1/18
Use geonet permanent station data for noise spectra analysis
in accordance with McNamara2004 using the obspy PPSD package
"""
import os
import sys
sys.path.append('./..')
import time
import glob
import argparse
from obspy import read, read_inventory
from obspy.clients.fdsn import Client
from obspy.signal import PPSD
from obspy import UTCDateTime
from obspy.imaging.cm import pqlx
from getdata import geonet_internal, fdsn_download, vog

# ignore warnings
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

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

num_of_days = int((end-start)/(3600*24))
print("\nTimeframe set: {} through {}, {} days".format(
                                        start.date,end.date,num_of_days))

# initialize PPSD with first datafile
oneday = 24*3600
dayplusone = start + oneday
start_check = time.time()
st, inv = fdsn_download(station=station,
                        channel=channel,
                        start=start,
                        end=dayplusone,
                        response=True)

if not st:
    sys.exit()

tr = st.select(channel=channel)[0]
if decimateby != 0:
    SR = tr.stats.sampling_rate
    print("SR {}->{}".format(SR,SR/decimateby))
    st.decimate(decimateby)
ppsd = PPSD(tr.stats, metadata=inv)
ppsd.add(st)
year_start = tr.stats.starttime.year
jday_start = tr.stats.starttime.julday
end_check = time.time()
print("complete ({}s)".format(round(end_check-start_check,2)))

# loop over rest of datafiles and add to ppsd
today = dayplusone
i=1
while today < end:
    try:
        print("{}/{}".format(i,num_of_days),end=' ')
        i+=1

        start_check = time.time()
        today = dayplusone
        dayplusone += oneday
        st, inv = fdsn_download(station=station,
                                channel=channel,
                                start=today,
                                end=dayplusone,
                                response=False)
        if decimateby != 0:
            st.decimate(decimateby)
        ppsd.add(st)
        end_check = time.time()
        print("complete ({}s)".format(round(end_check-start_check,2)))
    except Exception as e:
        print(e)
        pass

year_end = st[0].stats.starttime.year
jday_end = st[0].stats.starttime.julday

output_filename = '{sta}.{cha}.{yearS}.{day_start:0>3}-{yearE}.{day_end:0>3}'.format(
                                                        sta=station,
                                                        cha=channel,
                                                        yearS=year_start,
                                                        yearE=year_end,
                                                        day_start=jday_start,
                                                        day_end=jday_end)

npz_filepath = '../ppsd_arrays/{}.npz'.format(output_filename)
png_filepath = '../output_plots/ppsd_plots/{}.png'.format(output_filename)


# ignore overwrite protection
ppsd.save_npz(npz_filepath)
print("saved .npz file: ",npz_filepath)

ppsd.plot(filename=png_filepath,cmap=pqlx,show_mean=True)
print("saved .png file: ",png_filepath)
