"""
Module for getting data via internal pathname grabbing (@GNS) or with
FDSN webservices via obspy (@VIC)
"""
import os
import glob
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

def pathnames():
    """
    simplify pathname calling between computers, call function to return the
    correct pathname instead of setting them in every script
    """
    cwd = os.getcwd()
    if cwd == "/Users/chowbr/Documents/subduction/spectral":
        where = "GNS"
    elif cwd == "/seis/prj/fwi/bchow/spectral":
        where = "VIC"
    else:
        where = "OTHER"
    path_dictionary = {"spectral":os.path.join(cwd,''),
                        "rdf":os.path.join(cwd,'..','RDF_Array',''),
                        "plots":os.path.join(cwd,'output_plots',''),
                        "ppsd":os.path.join(cwd,'ppsd_arrays',''),
                        "where": where
                        }

    return path_dictionary

def geonet_internal(station,channel,start,end=False,response=True):
    """
    returns a list of pathnames for GEONET archives on GNS internal system.
    If response == True, also returns path for response.
    If end not specified, returns a list of length 1 for day requested.
    Wildcards possible, however inventory returns a full list of instruments at
    a certain location, rather than a paired down list.

    :type station: str
    :param station: station name i.e. GKBS (case-insensitive)
    :type channel: str
    :param channel: channel of interest, i.e. BHN, HHE (case-insensitive)
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
    station = station.upper()
    channel = channel.upper()
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
        print("Station choice does not exist")
        return False, False
    if not glob.glob(mseed_GNApath):
        print("Channel choice does not exist")
        return False, False
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
    # print("\n","="*25)
    # print("++ {nfiles} files for {days_between} days requested".format(
                                                    # nfiles=len(mseed_files),
                                                    # days_between=days_between))

    # response information; mseed sets naming parameters
    NET,STA,LOC,CHA,SUFX,YEAR,JDAY = os.path.basename(
                                                mseed_files[0]).split('.')
    if response:
        response_filename = "RESP.{net}.{sta}.{loc}.*".format(net=NET,
                                                            sta=STA,
                                                            loc=LOC)
        response_filepath = glob.glob(
                                os.path.join(resp_GNApath,response_filename))
        # print("++ response filepath")
    else:
        response_filepath = None

    # print(" ","="*25,"\n")

    return mseed_files, response_filepath

def fdsn_download(station,channel,start,end=False,response=False,
                                                    client="GEONET",cushion=0):
    """Download data via FDSN client for given station, channel and start
    and end times. Can output response as well. Return stream and response.
    :type station: str
    :param station: station name i.e. GKBS (case-insensitive)
    :type channel: str
    :param channel: channel of interest, i.e. BHN, HHE (case-insensitive)
    :type start: str
    :param start: starttime for data request
    :type end: str
    :param end: endtime for data request
    :type response: bool
    :param response: whether or not to return the response information, if not,
                    just returns bool False
    :type client: str
    :param client: client that fdsn searches for data
    :type cushion: int
    :param cushion: padding on start and end time
    """
    station = station.upper()
    folders = 1 # corresponds to number of channels
    if channel[-1] == "*":
        folders = 3
    start = UTCDateTime(start) - cushion
    if end:
        end = UTCDateTime(end) + cushion
    else:
        end = start + 3600*24

    # set instrument id from arguments
    instrument_id = 'NZ.{}.*.{}'.format(station,channel)
    # print("++ Requesting data for instrument: {}".format(instrument_id))
    net, sta, loc, cha = instrument_id.split('.')

    # initiate downloading of data
    c = Client(client)
    err_num = 0
    try:
        st = c.get_waveforms(network=net,
                            station=sta,
                            location=loc,
                            channel=cha,
                            starttime=start,
                            endtime=end)

    except Exception as e:
        print(e)
        err_num += 1
        if err_num > 5:
            print("errored out")
            a=1/0
        pass

    if response:
        # fetch response file
        # print("++ Requesting response information")
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

def event_stream(station,channel,event_id=False,pad=False):
    """Given a GEONET event ID, return raw waveform streams and response files
    Waveforms can be from RDF or GEONET permanent stations, chooses
    correct downloading format based on requests. If no event_id given, full
    catalog downloaded for times set in the function.

    :type choice: str
    :param choice: GEONET or RDF (temporary) station data
    :type station: str
    :param station: station name i.e. GKBS or RD01 (case-insensitive)
    :type channel: str
    :param channel: channel of interest, i.e. BHN, HHE (case-insensitive)
    :type event_id: str
    :param event_id: GEONET event id for earthquakes, if default, function will
                    simply collect a catalog of events for a given time periods
    :type pad: int
    :param pad: time in seconds to pad beginning of data, should be positive
    """
    # parse arguments
    station_id = station.upper() # i.e. RD06 or PXZ
    choice = "geonet"
    if station_id[:2] == "RD":
        choice = "rdf"
    channel = channel.upper()

    # get event information
    c = Client("GEONET")
    if event_id:
        cat = c.get_events(eventid=event_id)
    else:
        t_start = UTCDateTime("2017-320")
        t_end = UTCDateTime("2018-001")
        cat = c.get_events(starttime=t_start,
                            endtime=t_end,
                            minmagnitude=4,
                            maxmagnitude=6,
                            minlatitude=-50,
                            maxlatitude=-35,
                            minlongitude=165,
                            maxlongitude=180,
                            orderby="magnitude")

    for event in cat:
        origin = event.origins[0].time
        if pad:
            origin = origin - abs(pad)

        # temporary station filepaths
        if choice == "rdf":
            d1 = pathnames()['rdf'] + "July2017_Sep2017/DATA_ALL/"
            d2 = pathnames()['rdf'] + "Sep2017_Nov2017/DATA_ALL/"
            d3 = pathnames()['rdf'] + "Nov2017_Jan2018/DATA_ALL/"
            date_search = "{year}.{jday}".format(
                                            year=origin.year,jday=origin.julday)
            path_vert,path_north,path_east = [],[],[]
            for search in [d1,d2,d3]:
                path_vert += glob.glob(search+"{sta}.{date}*HHZ*".format(
                                                            sta=station_id,
                                                            date=date_search))
                path_north += glob.glob(search+"{sta}.{date}*HHN*".format(
                                                            sta=station_id,
                                                            date=date_search))
                path_east += glob.glob(search+"{sta}.{date}*HHE*".format(
                                                            sta=station_id,
                                                            date=date_search))

            st = read(path_vert[0]) + read(path_north[0]) + read(path_east[0])
            resp_filepath = pathnames()['rdf'] + "DATALESS.RDF.XX"
            inv = read_inventory(resp_filepath)


        # grab geonet data either internally or via fdsn
        elif choice == "geonet":
            if pathnames()["where"] == "GNS":
                path_vert, resp_vert = geonet_internal(station=station_id,
                                            channel= channel[:2] + 'Z',
                                            start = origin,
                                            response = True)
                path_north, resp_north = geonet_internal(station=station_id,
                                            channel= channel[:2] + 'N',
                                            start = origin,
                                            response = True)
                path_east, resp_east = geonet_internal(station=station_id,
                                            channel= channel[:2] + 'E',
                                            start = origin,
                                            response = True)
                inv = read_inventory(resp_vert)
                inv += read_inventory(resp_north)
                inv += read_inventory(resp_east)
                st = (read(path_vert[0]) + read(path_north[0]) +
                                                            read(path_east[0]))

            # if not working on GNS computer
            else:
                st, inv = fdsn_download(station=station_id,
                                        channel = channel[:2] + '*',
                                        start = origin,
                                        response = True)

        return st, inv, cat
