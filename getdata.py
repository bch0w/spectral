"""Downloads and saves data from FDSN webservice for instruments on FDSN,
over a given time period. Saves to folder under the current working directory in
a schema similar to observatory practices
i.e. ./NZ/'Station'/'Channel'/'data' -> ./NZ/TUDS/BHN/*.mseed

Example calls:
python getdata_fdsn.py --station GKBS --channel BN* --start 2016-11-13 \
--end 2016-11-14 --response True

"""
import os
import sys
import glob
import argparse
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

def pathnames(choice):
    """
    simplify pathname calling between Vic and GNS, call function to return the
    correct pathname instead of setting them in every script
    """
    choice = choice.upper()
    if choice == 'VIC':
        basepath = '/Users/chowbr/Documents/subduction/'
        path_dictionary = {"spectral":basepath + 'spectral/',
                            "rdf":basepath + 'RDF_Array/',
                            "plots":basepath + 'spectral/output_plots/',
                            "ppsd":basepath + 'spectral/ppsd_arrays/'}
    elif choice == 'GNS':
        basepath = '/seis/prj/fwi/bchow/'
        path_dictionary = {"spectral":basepath + 'spectral/',
                            "rdf_y":'/seis/prj/fwi/yoshi/RDF_Array/',
                            "rdf_b":basepath + 'RDF_Array/',
                            "plots":basepath + 'spectral/output_plots/',
                            "ppsd":basepath + 'ppsd_arrays/'}

    return path_dictionary

def geonet_internal(station,comp,start,end=False,response=True):
    """
    returns a list of pathnames for GEONET archives on GNS internal system.
    If response == True, also returns path for response.
    If end not specified, returns a list of length 1 for day requested.
    Currently only works for single component requests, no wildcards possible.

    :type station: str
    :param station: station name i.e. GKBS (case-insensitive)
    :type comp: str
    :param comp: component of interest, N/E/Z (case-insensitive)
    :type start: str
    :param start: starttime for data request
    :type end: str
    :param end: endtime for data request
    :type response: bool
    :param response: return response filepath
    :rtype mseed_files: list of str
    :return mseed_files: list of absolute filepaths to requested data
    :rtype response_filepath: str or None
    :return response_filepath: if requested, filepath of response
    """
    # channel naming convention based on instrument type
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

    # check if station exists
    check_sta_path = '/geonet/seismic/{year}/NZ/'.format(year=start.year)
    if not glob.glob(check_sta_path + station):
        sys.exit("Station choice does not exist")
    if not glob.glob(mseed_GNApath):
        sys.exit("Channel choice does not exist")

    # check if data spans more than one day
    if type(end) != bool:
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
        else:
            # first year
            first_year_files = glob.glob(mseed_GNApath + '*')
            first_year_files.sort()
            start_file_match = glob.glob(mseed_GNApath + "*{date}".format(
                                date=start.format_seed().replace(',','.')))[0]
            start_file_index = first_year_files.index(start_file_match)
            mseed_files = first_year_files[start_file_index:]
            # if years in middle
            if end.year - start.year > 1:
                for year_iter in range(start.year+1,end.year,1):
                    pth_iter = '/geonet/seismic/{year}/NZ/{sta}/{ch}.D/'.format(
                                                                year=year_iter,
                                                                sta=station,
                                                                ch=channel)
                    mseed_files += glob.glob(pth_iter + '*')
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

        days_between = int((end-start)/86400)
    # if data only needed for one day
    else:
        mseed_files = glob.glob(mseed_GNApath +'*{year}.{day:0>3}'.format(
                                                            year=start.year,
                                                            day=start.julday))
        days_between = 1

    mseed_files.sort()
    print("\n","="*25)
    print("++ {nfiles} files for {days_between} days requested".format(
                                                    nfiles=len(mseed_files),
                                                    days_between=days_between))

    # response information; mseed sets naming parameters
    NET,STA,LOC,CHA,SUFX,YEAR,JDAY = os.path.basename(
                                                mseed_files[0]).split('.')
    if response:
        response_filename = "RESP.{net}.{sta}.{loc}.{cha}".format(net=NET,
                                                                    sta=STA,
                                                                    loc=LOC,
                                                                    cha=CHA)
        response_filepath = os.path.join(resp_GNApath,response_filename)
        print("++ response filepath")
    else:
        response_filepath = None

    print(" ","="*25,"\n")

    return mseed_files, response_filepath

def fdsn_download(station,channel,start,end,response=False,client="GEONET",):
    """Download data via FDSN client for given station, channel and start
    and end times. Can output response as well. Return stream and response.

    """
    station = station.upper()
    folders = 1 # corresponds to number of channels
    if channel[-1] == "*":
        folders = 3
    start = UTCDateTime(start)
    end = UTCDateTime(end)

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
        # fetch waveforms
        try:
            st = c.get_waveforms(network=net,
                                station=sta,
                                location=loc,
                                channel=cha,
                                starttime=day_of,
                                endtime=day_of+sec_per_day)
            day_of += sec_per_day

        except Exception as e:
            print(e)
            err_num += 1
            if err_num > 5:
                sys.exit('Errored out')
            pass

    if response:
        # fetch response file
        print("++ Requesting response information")
        try:
            response = c.get_stations(network=net,
                                station=sta,
                                location='*',
                                channel='*',
                                starttime=start,
                                endtime=end,
                                level='response')
        except Exception as e:
            print(e)
            pass


    return st, response

if __name__ == "__main__":
    # user input arguments
    parser = argparse.ArgumentParser(description='Download data via FDSN client\
     GEONET, example call \n $ python getdata_fdsn.py --station TBAS --channel\
     BNZ --end 2015-02-01 --response True')
    parser.add_argument('--station', help='Instrument station, default = KNZ',
                        type=str,default='KNZ')
    parser.add_argument('--channel', help='Instrument channel, i.e. BN1/HHZ, \
                        if ??* downloads all 3 components, default = HHZ',
                        type=str, default='Z')
    parser.add_argument('--start', help='Starttime YYYY-MM-DD default = \
                        2015-01-01',type=str, default='2015-01-01')
    parser.add_argument('--end', help='Endtime, default = 2015-01-02',type=str,
                        default='2015-01-02')
    parser.add_argument('--response', help='Boolean, default = False',type=bool,
                        default=False)
    parser.add_argument('--client', help='FDSN client, default = GEONET',
                        type=str, default='GEONET')

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
                                    '*{year}.{day:0>3}*'.format(
                                                            year=day_of.year,
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
