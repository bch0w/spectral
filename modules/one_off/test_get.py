
import os
import sys
import glob
import argparse
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

station = 'gkbs'
comp = 'z'
start = '2015-01-01'
end = '2017-02-01'
# end = False
response = True

chan_dict = {'Z':'HH','S':'BN'}
station = station.upper()
comp = comp.upper()
channel = chan_dict[station[-1]] + comp # i.e. HHZ
start = UTCDateTime(start)

# filepaths direct to GeoNet Archives path
mseed_GNApath = '/geonet/seismic/{year}/NZ/{sta}/{ch}.D/'.format(
                                                            year=start.year,
                                                            sta=station,
                                                            ch=channel)
resp_GNApath = '/geonet/seed/RESPONSE/{sta}.NZ/'.format(sta=station)

# check if data spans more than one day
if end:
    end = UTCDateTime(end)
    # if data only spans one year
    if start.year == end.year:
        list_of_files = glob.glob(mseed_GNApath + '*')
        list_of_files.sort()
        # match start and end dates to filepaths
        start_file_match = glob.glob(mseed_GNApath + "*{date}".format(
                                date=start.format_seed().replace(',','.')))[0]
        end_file_match = glob.glob(mseed_GNApath + "*{date}".format(
                                date=end.format_seed().replace(',','.')))[0]
        # determine indices for start and end files
        start_file_index = list_of_files.index(start_file_match)
        end_file_index = list_of_files.index(end_file_match)
        # list of files for return
        mseed_files = list_of_files[start_file_index:end_file_index]

    # if data spans multiple years
    if start.year != end.year:
        # first year
        first_year_files = glob.glob(mseed_GNApath + '*')
        first_year_files.sort()
        start_file_match = glob.glob(mseed_GNApath + "*{date}".format(
                                 date=start.format_seed().replace(',','.')))[0]
        start_file_index = first_year_files.index(start_file_match)
        mseed_files = first_year_files[start_file_index:]
        # if years in middle
        if end.year - start.year > 1:
            for year_iter in range(start.year,end.year,1):
                GNApath_iter = '/geonet/seismic/{year}/NZ/{sta}/{ch}.D/'.format(
                                                                year=year_iter,
                                                                sta=station,
                                                                ch=channel)
                mseed_files += glob.glob(GNApath_iter + '*')
        # last year
        GNApath_last = '/geonet/seismic/{year}/NZ/{sta}/{ch}.D/'.format(
                                                                year=end.year,
                                                                sta=station,
                                                                ch=channel)
        last_year_files = glob.glob(GNApath_last + '*')
        last_year_files.sort()
        end_file_match = glob.glob(GNApath_last + "*{date}".format(
                                 date=end.format_seed().replace(',','.')))[0]
        end_file_index = last_year_files.index(end_file_match)
        mseed_files += last_year_files[:end_file_index]

# if data only needed for one day
else:
    mseed_files = glob.glob(mseed_GNApath +'*{year}.{day:0>3}'.format(
                                                        year=start.year,
                                                        day=start.julday))
mseed_files.sort()

# response information; mseed sets naming parameters
NET,STA,LOC,CHA,SUFX,YEAR,JDAY = os.path.basename(
                                            mseed_files[0]).split('.')
if response:
    response_filename = "RESP.{net}.{sta}.{loc}.{cha}".format(net=NET,
                                                                sta=STA,
                                                                loc=LOC,
                                                                cha=CHA)
    response_filepath = os.path.join(resp_GNApath,response_filename)
