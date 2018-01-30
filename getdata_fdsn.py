"""Downloads data from FDSN webservice for GEONET permanent stations,
over a given time period. Downloads data for full days, saves to folder under
the current working directory in a schema similar to observatory practices
i.e. ./NZ/'Station'/'Channel'/'data' -> ./NZ/TUDS/BHN/*.mseed

!!! If you're working internally at GNS, the data can be accessed under the
!!! directory /geonet/seismic/ and the response files can be accessed under the
!!! directory /geonet/seed/RESPONSE

Example calls:
python getdata_fdsn.py --station GKBS --channel BN* --start 2016-11-13 --end 2016-11-14

"""
import os
import sys
import glob
import argparse
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# user input arguments
parser = argparse.ArgumentParser(description='Download data via FDSN client \
GEONET, example call \n $ python getdata_fdsn.py --station TBAS --channel BNZ \
--end 2015-02-01 --response True')
parser.add_argument('--station', help='Instrument station, default = KNZ',type=str,
                        default='KNZ')
parser.add_argument('--channel', help='Instrument channel, i.e. BN1/HHZ, if ??*\
                    downloads all 3 components, default = HHZ',type=str,
                    default='Z')
parser.add_argument('--start', help='Starttime YYYY-MM-DD default = 2015-01-01',
                    type=str, default='2015-01-01')
parser.add_argument('--end', help='Endtime, default = 2015-01-02',type=str,
                    default='2015-01-02')
parser.add_argument('--response', help='Boolean, default = False',type=bool,
                    default=False)
parser.add_argument('--client', help='FDSN client, default = GEONET',type=str,
                    default='GEONET')

# parse arguments
arg = parser.parse_args()
station = arg.station.upper()
channel = arg.channel
folders = 1 # corresponds to number of channels
if channel[-1] == "*":
    folders = 3
start = UTCDateTime(arg.start)
end = UTCDateTime(arg.end)
response = arg.response
client = arg.client

# set instrument id from arguments
instrument_id = 'NZ.{}.*.{}'.format(station,channel)
print("++ Requesting data for intrument: {}".format(instrument_id))
net, sta, loc, cha = instrument_id.split('.')

# split time into days
sec_per_day = 3600*24
time_delta_days = int((end-start)/(sec_per_day))

# initiate downloading of data
c = Client(client)
day_of = start
err_num = 0
while day_of <= end:
    # check if data exists, if so, iterate on day
    data_check = os.path.join(net,sta,'{}'.format(channel),
                                '*{year}.{day:0>3}*'.format(year=day_of.year,
                                                            day=day_of.julday))
    if len(glob.glob(data_check)) == folders:
        print(data_check,' already exists')
        day_of += sec_per_day
        continue

    # check if path for station already exists already exists
    filepath_data = os.path.join(net,sta)
    if not os.path.exists(filepath_data):
         os.makedirs(filepath_data)

    # fetch waveforms
    try:
        st = c.get_waveforms(network=net,
                            station=sta,
                            location=loc,
                            channel=cha,
                            starttime=day_of,
                            endtime=day_of + sec_per_day)

        # address all traces in stream
        for tr in st:
            print(tr)
            # # if previous day captured, skip over
            # if tr.stats.starttime.julday != day_of.julday:
            #     print("\t skipped")
            #     continue
            stats = tr.stats
            # check path, if not, then make one
            channel_path = os.path.join(filepath_data,stats.channel)
            if not os.path.exists(channel_path):
                os.makedirs(channel_path)

            # create file naming and folder naming schema
            output_data_name = '{id}.{year}.{jday:0>3}.mseed'.format(
                            id=tr.id,
                            year=day_of.year,
                            jday=day_of.julday)
            full_data_path = os.path.join(channel_path,output_data_name)

            tr.write(full_data_path,format="MSEED")
            print('\t',output_data_name," written")

        # iterate on day
        day_of += sec_per_day

    except Exception as e:
        print(e)
        err_num += 1
        if err_num > 5:
            sys.exit('Errored out')
        pass

if response:
    print("++ Requesting response information")
    output_inv = '{net}_{sta}.xml'.format(net=net,sta=sta)
    filepath_inv = os.path.join(net,sta,output_inv)

    # check if inventory already fetched
    if len(glob.glob(filepath_inv)) == 1:
        sys.exit('Inventory exists')
    try:
        inv = c.get_stations(network=net,
                            station=sta,
                            location='*',
                            channel='*',
                            starttime=start,
                            endtime=end,
                            level='response')
        inv.write(filepath_inv,format="STATIONXML")
        print(output_inv,' written')
    except Exception as e:
        print(e)
        pass
