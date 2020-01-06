"""
Download data from FDSN based on GCMT event ids, write data to specified format
Data will be saved per channel. Simple preprocessing available

This requires ObsPy (obspy.org) and Python3
"""
import os
import numpy as np
from obspy import UTCDateTime, read_events
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth


def query_gcmt(event_id, return_idx=0):
    """
    Query the GCMT online web service for events matching a GCMT event id
    If no matching event ID found, will search based on time criteria set in
    the GCMT event id. If nothing found returns None

    Note: assumes the GCMT event id is defined by origintime 
    This is a rewritten version of: 
        pyatoa.utils.gather.grab_auxiliaries.grab_gcmt_moment_tensor

    :type event_id: str
    :param event_id: GCMT event id, e.g. "C200902170330A"
    :type return_idx: int
    :param return_idx: special use case; if no event by event id is found, but
        multiple events are found using origin time, this index determines
        which event in the catalog to choose
    :rtype: obspy.core.event.event.Event
    :return: obspy Event object which contains relevant information for event
    """
    from urllib.error import HTTPError

    gcmt_monthly= ("https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/"
                   "NEW_MONTHLY/{year_long}/{month}{year_short}.ndk")
    gcmt_quick = ("https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/"
                  "NEW_QUICK/gcmt.ndk")
  
    # Get origin time information for query
    # GCMT event ids can come prepended with a C or appended with A, strip off
    origin_str = event_id.replace("C", "")
    origin_str = event_id.replace("A", "")
    origin_time = UTCDateTime(origin_str)

    # Date strings for HTML parsing
    month = origin_time.strftime('%b').lower()
    year_short = origin_time.strftime('%y')
    year_long = origin_time.strftime('%Y')

    # Query the GCMT website using Obspy
    try:
        cat = read_events(gcmt_monthly.format(year_long=year_long, 
                                              year_short=year_short,
                                              month=month))
    except HTTPError:
        cat = read_events(gcmt_quick)

    # Try to search the catalog for the given event id
    for event in cat:
        if event_id in event.resource_id.id:
            return event
   
    # If event id query failed, try to filter the catalog by origintime +/-100s    
    cat_filt = cat.filter("time > {}".format(str(origintime - 100)),
                          "time < {}".format(str(origintime + 100))
                          )
    # The filter returns Catalogs, filter based on length
    if bool(cat_filt):
        if len(cat_filt) == 1:
            return cat_filt
            import ipdb;ipdb.set_trace()
        else:
            print("multiple events, returning index {}".format(return_idx))
            return cat_filt[return_idx]
    else:
        print("No events found for event id {}".format(event_id))
        return


if __name__ == "__main__":
    # ==========================  SET PARAMETERS BELOW =========================

    # Station
    stations_file = "./stations_gsn.txt"  # path to station list
    channel = "BH*"  # channels to gather, e.g. "BHZ", "*", "?H?", wildcards ok

    # Event
    id_or_time = "id"  # Picks between two lists below, options: 'id' or 'time'
    event_ids = ["200902170330A", "201712080209A"]  # Pick between GCMT event id
    origin_times = []  # OR origin times e.g. "2009-02-17T03:30:58.6"
    seismogram_length_seconds = 60 * 60  # length of seismo after origin

    # Processing
    process = True  # can skip all processing if you just want raw seismograms
    remove_response = True  # if you want to remove the instrument response
    min_period = 10  # None for highpass. Bandpass if max_period also set
    max_period = None # None for lowpass. Bandpass if min_period also set

    # Output
    output_directory = "./"  # directory to save waveform data to
    output_format = "SAC"  # 'MSEED', 'SAC', 'SEGY', 'SU' etc., see ObsPy 
    save_raw_seismograms = False  # if process == True, also save raw seismo
    
    # ==========================  SET PARAMETERS ABOVE =========================

    failed = 0
    c = Client("IRIS")
    # If event id's chosen, query GCMT to get starttimes
    events, event_distances = [], []
    if id_or_time == "id":
        origin_times = []
        for eid in event_ids:
            event = query_gcmt(eid)
            events.append(event)
            event_distances.append(eid + "_dists.txt")
            origin_times.append(event.preferred_origin().time)

    if event_distances:
        for edid in event_distances:
            if not os.path.exists(edid):
                with open(edid, 'w') as f:
                    f.write("NET,STA,DIST(km),DIST(deg),AZ(deg),BAz(deg)\n")

    # Use starttimes to define events
    for i, starttime in enumerate(origin_times):
        starttime = UTCDateTime(starttime)
        endtime = starttime + seismogram_length_seconds 
        print(starttime)
        
        stations = np.loadtxt(stations_file, dtype="str")
        for station_list in stations:
            network = station_list[1]
            station = station_list[0]

            # Try calculate the azimuth and distance
            if event_distances:
                dist_m, az, baz = gps2dist_azimuth(
                                    lat1=events[i].preferred_origin().latitude,
                                    lon1=events[i].preferred_origin().longitude,
                                    lat2=float(station_list[2]),
                                    lon2=float(station_list[3])
                                    )
                dist_deg = (dist_m * 1E-3) / 111.11
                with open(event_distances[i], "a") as f:
                    f.write(
                     "{net},{sta},{dst_km:.2e},{dst_deg:.2f},{az:.2f},{baz:6.2f}\n".format(
                                                     net=network, sta=station,
                                                     dst_km = dist_m * 1E-3,
                                                     dst_deg=dist_deg, az=az, 
                                                     baz=baz)
                                                     )
            # Set the filename for saving the waveforms
            fid_out = os.path.join(output_directory, "{time}_{net}_{sta}".format(
                                   net=network, sta=station, 
                                   time=str(starttime).split('T')[0])
                                   )

            # Get the waveform data from FDSN
            print(network, station, end="... ")
            try:
                st = c.get_waveforms(network=network, station=station, 
                                     location="*", channel=channel, 
                                     starttime=starttime, endtime=endtime,
                                     attach_response=remove_response)

                # Write the raw data
                if save_raw_seismograms:
                    for tr in st:

                        fid_out_raw = fid_out + "_{cha}_raw.{fmt}".format(
                                                     cha=tr.stats.channel, 
                                                     fmt=output_format.lower())
                        tr.write(fid_out_raw, format=output_format)
                if process:
                    # Remove response
                    if remove_response:
                        st.remove_response()

                    # Preprocessing functionality
                    st.detrend("linear")  # detrend the data, remove long period
                    st.detrend("demean")  # demean to set data to 0
                    st.taper(max_percentage = 0.05)  # taper the ends 

                    # Filter the data 
                    if min_period and max_period:
                        st.filter("bandpass", freqmin=1/max_period, 
                                  freqmax=1/min_period)
                    elif min_period and not max_period:
                        st.filter("lowpass",  freq=1/min_period)
                    elif max_period and not min_period:
                        st.filter("highpass", freq=1/max_period)

                    st.detrend("linear")  # detrend and taper again to remove
                    st.detrend("demean")  # introduced spurious signals
                    st.taper(max_percentage = 0.05)

                    # Save the data based on the User-defined file format
                    for tr in st:

                        fid_out_proc = fid_out + "_{cha}_proc.{fmt}".format(
                                                     cha=tr.stats.channel,
                                                     fmt=output_format.lower())
                        tr.write(fid_out_proc, format=output_format)
                print("gathered")
            except Exception as e:
                failed += 1
                print("failed")
                print(e)
                continue





