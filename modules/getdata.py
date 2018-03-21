"""
Module for getting data via internal pathname grabbing (@GNS) or with
FDSN webservices via obspy (@VIC)
"""
import os
import glob
from obspy import read, read_events, read_inventory, UTCDateTime
from obspy.clients.fdsn import Client

# ============================== INTERNAL SCRIPTS =============================
def pathnames():
    """
    simplify pathname calling between computers, call function to return the
    correct pathname instead of setting them in every script
    """
    cwd = os.getcwd()
    base0 = cwd.split('/')[0]
    base1 = cwd.split('/')[1]
    base2 = cwd.split('/')[2]
    basecheck = os.path.join(base0,base1,base2)
    # if cwd == "/Users/chowbr/Documents/subduction/spectral":
    #     where = "VIC"
    # elif cwd == "/seis/prj/fwi/bchow/spectral":
    #     where = "GNS"
    # else:
    #     where = "OTHER"
    if basecheck == "seis/prj":
        where = "VIC"
        base = "/seis/prj/fwi/bchow/spectral"
    elif basecheck == "Users/chowbr":
        where = "VIC"
        base = "/Users/chowbr/Documents/subduction/spectral"

    common = os.path.join(base,'common')
    path_dictionary = {"spectral":os.path.join(base,'spectral',''),
            "kupe":os.path.join(base,'kupe',''),
            "rdf":os.path.join(common,'DATA','RDF_Array',''),
            "plots":os.path.join(common,'OUTPUT_PLOTS',''),
            "spectralplots":os.path.join(common,'OUTPUT_PLOTS','spectral',''),
            "kupeplots":os.path.join(common,'OUTPUT_PLOTS','kupe',''),
            "ppsd":os.path.join(common,'DATA','ppsd_arrays',''),
            "data":os.path.join(common,'DATA',''),
            "kupedata":os.path.join(common,'DATA','KUPEDATA',''),
            "hobitss":os.path.join(common,'DATA','hobitss_mseeds',''),
            "mseed":os.path.join(common,'DATA','MSEED',''),
            "syns":os.path.join(common,'DATA','MSEED','SYN',''),
            "where": where
                    }

    return path_dictionary

# ============================== WAVEFORM DOWNLOAD =============================

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
        response_filename = "RESP.{net}.{sta}.{loc}.{cha}".format(net=NET,
                                                            sta=STA,
                                                            loc=LOC,
                                                            cha=CHA)
        response_filepath = glob.glob(
                                os.path.join(resp_GNApath,response_filename))[0]
        # print("++ response filepath")
    else:
        response_filepath = None

    # print(" ","="*25,"\n")

    return mseed_files, response_filepath

def fdsn_download(code,event_id=False,response=False,
                                        client="GEONET",cushion=500):
    """Download data via FDSN client for given station, channel and start
    and end times. Can output response as well. Return stream and response.
    :type code: str
    :param code: instrument code in the format NN.SSSS.LL.CCC i.e. NZ.PUZ..HHZ
    :type start: str
    :param start: starttime for data request
    :type end: str
    :param end: endtime for data request, default 2 hours after start + cushion
    :type response: bool
    :param response: whether or not to return the response information, if not,
                    just returns bool False
    :type client: str
    :param client: client that fdsn searches for data
    :type cushion: int
    :param cushion: padding on start and end time, default 500s
    """
    # set timing
    if event_id:
        cat = get_quakeml(event_id)
        event = cat[0]
        starttime = event.origins[0].time - cushion
        endtime = starttime + 3600*2 + cushion
    elif start:
        starttime = UTCDateTime(start) - cushion
        if end:
            endtime = UTCDateTime(end) + cushion
        else:
            endtime = start + 3600*2 + cushion

    # set instrument id from arguments
    net,sta,loc,cha = code.split('.')

    # special sets for hobitss data
    if net == 'YH':
        # skip any event that falls outside the hobitss timeframe
        global_start = UTCDateTime("2014-05-13T12:49:00.000000Z")
        global_end = UTCDateTime("2015-06-20T13:55:00.000000Z")
        if (starttime <= global_start) or (endtime >= global_end):
            print(eventid,"not within Hobitss timeframe")
            return None, None
        client = "IRIS"

    # search for data internal
    data_path = pathnames()['mseed'] + '{n}/{e}_{s}.mseed'.format(n=net,
                                                                 e=event_id,
                                                                 s=sta)
    try:
        st = read(data_path)
    except Exception as e:
        print("[getdata.fdsn_download] internal mseed not found,"
                                                            " querying fdsn...")

        # initiate downloading of data
        c = Client(client)
        err_num = 0
        try:
            st = c.get_waveforms(network=net,
                                station=sta,
                                location=loc,
                                channel=cha,
                                starttime=starttime,
                                endtime=endtime)
            st.write(data_path,format="MSEED")
            print('\t++',data_path)
        except Exception as e:
            print(e)

    if response:
        # check existing stationxml files
        response_dict = {"NZ":"NZ_NI_BB_seismometers.xml",
                         "YH":"YH_LOBS.xml"}

        resp_path = (pathnames()['data'] +
                                    'STATIONXML/{}'.format(response_dict[net]))
        inv_full = read_inventory(resp_path)
        inv_check = inv_full.select(station=sta)

        if len(inv_check) > 0:
            inv = inv_check
        else:
            print("[getdata.fdsn_download] internal resp not found,"
                                                            " querying fdsn...")
            try:
                c = Client(client)
                inv_get = c.get_stations(network=net,
                                        station=sta,
                                        location=loc,
                                        channel=cha,
                                        starttime=starttime,
                                        endtime=endtime,
                                        level='response')
                inv_full += inv_get[0]
                # ============================================================
                # !!! writes into a new network, even if network is the same
                # !!! okay if you use the inv.select command, but
                # !!! aesthetically annoying
                inv_full.write(resp_path,format="STATIONXML")
                # ============================================================
                inv = inv_get
                print('\t++',resp_path)
            except Exception as e:
                print(e)
                return st, None

        return st, inv


    else:
        return st, None


def event_stream(code,event_id,client="GEONET",startpad=False,endpad=False):
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
    :type startpad: int
    :param startpad: time in seconds to pad beginning of data, positive
    :type endpad: int
    :param endpad: time in seconds from starttime to end cutoff
    """
    net,sta,loc,cha = code.split('.')
    cat = get_quakeml(event_id)

    # determine if RDF data or other
    choice = "NOTRDF"
    if sta[:2] == "RD":
        choice = "RDF"

    for event in cat:
        # set start and end times
        origin = event.origins[0].time
        if startpad:
            start = origin - startpad
        else:
            start = origin

        if endpad:
            end = start + endpad
        else:
            end = start + 1000

        # grab geonet data either internally or via fdsn
        if choice == "NOTRDF":
            if pathnames()["where"] == "GNS":
                st_list,resp_list = [],[]
                for comp in ['N','E','Z']:
                    path_st, path_resp = geonet_internal(station=sta,
                                            channel=cha[:2] + comp,
                                            start=origin,
                                            response=True)
                    st_list += path_st
                    resp_list.append(path_resp)

                inv = (read_inventory(resp_list[0]) +
                        read_inventory(resp_list[1]) +
                        read_inventory(resp_list[2])
                        )
                st = (read(st_list[0]) +
                        read(st_list[1]) +
                        read(st_list[2])
                        )

            # if not working on GNS computer
            else:
                st, inv = fdsn_download(code=code,
                                        event_id=event_id,
                                        response=True)

            # trim regardless of get-method
            st.trim(starttime=start,endtime=end)

        # temporary station filepaths
        if choice == "RDF":
            st, inv = get_rdf_data(origin)

        return st, inv, cat

def get_rdf_data(origin):
    """get data from RDF array internally UNTESTED, needs work probably
    """
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

    return st, inv


# ============================== EVENT INFO DOWNLOAD ==========================

def get_quakeml(event_id):
    """get quakeml file for event information, save internally if not already
    present in filesystem
    """
    event_path = pathnames()['data'] + 'QUAKEML/{}.xml'.format(event_id)
    try:
        cat = read_events(event_path)
    except Exception as e:
        print("[getdata.get_quakeml] internal quakeml not found,"
                                                        " querying fdsn... ")
        event_client = Client("GEONET")
        cat = event_client.get_events(eventid=event_id)
        cat.write(event_path,format="QUAKEML")
        print('\t++',event_path)

    return cat

def get_moment_tensor(event_id):
    """gets moment tensor as array from geonet CSV file"""
    import csv
    csvfile = pathnames()['data'] + "GeoNet_CMT_solutions.csv"
    with open(csvfile,'r') as f:
        reader = csv.reader(f)
        for i,row in enumerate(reader):
            if i == 0:
                tags = row
            if event_id == row[0]:
                values = []
                for t,v in zip(tags,row):
                    if (t == "Date") or (t == "PublicID"):
                        values.append(v)
                    else:
                        values.append(float(v))

                MT = dict(zip(tags,values))
                return MT

def get_GCMT_solution(event_id):
    """retrieve GCMT solutions from the .ndk files for comparison against
    converted MT from Ristau. Returns an obspy event object.
    """
    import os
    from obspy import read_events, UTCDateTime
    from getdata import get_moment_tensor

    month_dict={4:"apr",12:"dec",1:"jan",6:"jun",5:"may",10:"oct",
                8:"aug",2:"feb",7:"jul",3:"mar",11:"nov",9:"sep"}
    MT = get_moment_tensor(event_id=event_id)
    mw = MT["Mw"]
    date = UTCDateTime(MT["Date"])
    year = str(date.year)
    month = date.month
    fid = "{m}{y}.ndk".format(m=month_dict[month],y=year[2:])
    filepath = os.path.join(pathnames()['data'],"GCMT",year,fid)
    cat = read_events(filepath)
    cat_filt = cat.filter("time > {}".format(str(date-300)),
                          "time < {}".format(str(date+300)),
                          "magnitude >= {}".format(mw-.5),
                          "magnitude <= {}".format(mw+.5)
                          )
    if len(cat_filt) == 0:
        print("No events found")
        return
    elif len(cat_filt) > 1:
        print("{} events found, choose from list:".format(len(cat_filt)))
        print(cat_filt)
        choice = input("Event number: ")
        event = cat_filt[choice]
        return event
    else:
        event = cat_filt[0]
        return event

# ============================== COPY-PASTE SCRIPTS ============================

def get_those_stations():
    """misc station getter to be copy-pasted"""
    from obspy.clients.fdsn import Client
    c = Client("GEONET")
    north_island = [-42,-34,173,180]
    north_island_zoom = [-40,-37,176,178.5]
    new_zealand = [-50,-35,165,180]
    lat_lon = new_zealand
    inv = c.get_stations(network='NZ',
                        station='*',
                        channel='HH*',
                        minlatitude=lat_lon[0],
                        maxlatitude=lat_lon[1],
                        minlongitude=lat_lon[2],
                        maxlongitude=lat_lon[3],
                        level="station")
    c = Client("IRIS")
    inv = c.get_stations(network='YH',
                        station="LOBS*",
                        location='',
                        level="station")
    inv += c.get_stations(network='YH',
                        station="EBS*",
                        location='',
                        level="station")

    inv.write('STATION_PLACEHOLDER.xml',format='STATIONXML')

def get_those_events():
    """misc event getter to be copy-pasted"""
    from obspy.clients.fdsn import Client
    c = Client("GEONET")
    # for full catalog
    cat = c.get_events(starttime=t_start,
                        endtime=t_end,
                        minmagnitude=4,
                        maxmagnitude=6,
                        minlatitude=-50,
                        maxlatitude=-35,
                        minlongitude=165,
                        maxlongitude=180,
                        orderby="magnitude")
    # for single event
    from obspy.clients.fdsn import Client
    c = Client("GEONET")
    event_id = "2015p822263"
    cat = c.get_events(eventid=event_id)
    # cat.write("EVENT_PLACEHOLDER.xml",format="QUAKEML")

def station_and_event():
    """copy paste script for plotting hobitss stations and event as beachball
    """
    import matplotlib.pyplot as plt
    from getdata import get_moment_tensor
    from obspy import read_inventory
    from obspy.clients.fdsn import Client
    from obspy.imaging.beachball import beach

    eventid="2014p864702"
    c = Client("GEONET")
    inv = read_inventory('./datafiles/hobitss_stations.xml')

    MT = get_moment_tensor(eventid)
    fig = inv.plot(continental_fill_color='w',
                    water_fill_color='w',
                    show=False)

    eventx,eventy = fig.bmap(MT['Longitude'],MT['Latitude'])
    FM = [MT['strike1'],MT['dip1'],MT['rake1']]
    b = beach(FM,xy=(eventx,eventy),width=3E4,linewidth=1,facecolor='r')
    b.set_zorder(10)
    ax = plt.gca()
    ax.add_collection(b)
    plt.annotate("M{}".format(MT['ML']),
                    xy=(eventx,eventy),
                    xytext=(eventx,eventy),
                    fontsize=12,
                    zorder=200,
                    weight='bold')
    plt.show()

if __name__ == "__main__":
    print('what?')
