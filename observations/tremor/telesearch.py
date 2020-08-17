"""function to find teleseismic events from GCMT catalog given a specific 
Julian day. Returns relevant information for FATHOM stations such as distance to
event, expected arrivals of surface waves etc.
"""
import os
import sys
import pytz
from datetime import timedelta
from obspy import UTCDateTime, read_events
from obspy.geodetics import gps2dist_azimuth

sys.path.append("../modules")
from getdata import pathnames


def search_GCMT(date):
    """lookup function for finding events in GCMT catalogs stored internally
    :type date: str format accepted by UTCDateTime
    """
    date = UTCDateTime(date)
    day_before = date - (60*60*24)
    two_days_before = date - (2*60*60*24)
    day_after = date + (60*60*24)
    
    month_dict={4:"apr",12:"dec",1:"jan",6:"jun",5:"may",10:"oct",
                8:"aug",2:"feb",7:"jul",3:"mar",11:"nov",9:"sep"}
    year = str(date.year)
    month = date.month
    
    fid = "{m}{y}.ndk".format(m=month_dict[month],y=year[2:])
    filepath = os.path.join(pathnames()['data'],"GCMT",year,fid)

    # files can also be read directly from GCMT website
    gcmt_standard_url = ("https://www.ldeo.columbia.edu/~gcmt/projects/CMT/"
                            "catalog/NEW_MONTHLY/{y}/"
                            "{m}{ys}.ndk".format(y=year,
                                                 m=month_dict[month],
                                                 ys=year[2:]))
    gcmt_quick_url = ("http://www.ldeo.columbia.edu/~gcmt/projects/CMT/"
                      "catalog/NEW_QUICK/qcmt.ndk")

    try:
        cat = read_events(filepath)
    except FileNotFoundError:
        try:
            print("[getdata.get_GCMT_solution] internal .ndk file not found, "
                  "searching for GCMT standard url")
            cat = read_events(gcmt_standard_url)
        except Exception as e:
            print("[getdata.get_GCMT_solution] standard url not found, "
                  "searching GCMT quick solutions")
            cat = read_events(gcmt_quick_url)
    cat_filt = cat.filter("time > {}".format(str(two_days_before)),
                          "time < {}".format(str(day_after)),
                          "magnitude >= {}".format(7.5)
                          )
    if len(cat_filt) == 0:
        print("[getdata.get_GCMT_solution] No events found")
        return None
    elif len(cat_filt) > 1:
        print("[getdata.get_GCMT_solution]"
                    " {} events found, choose from list:".format(len(cat_filt)))
        print(MT)
        print(cat_filt)
        choice = int(input("Event number (index from 0): "))
        event = cat_filt[choice]
        return event
    else:
        event = cat_filt[0]
        return event
        
def fathom_event(event):
    """print out parameters related to event and fathom
    """
    fthm_lat,fthm_lon = [-40.0643,176.441]
    evnt_lat,evnt_lon = event.origins[0].latitude,event.origins[0].longitude
    G2DA = gps2dist_azimuth(lat1=evnt_lat,lon1=evnt_lon,
                            lat2=fthm_lat,lon2=fthm_lon)
    great_circ_dist_m,Az,BAz = G2DA
    great_circ_dist_km = great_circ_dist_m * 1E-3
    approximate_gcd_deg = great_circ_dist_km/111.11
    
    # convert UTCtime to local-time
    startUTC = event.origins[0].time.datetime.replace(tzinfo=pytz.utc)
    startLOC = startUTC.astimezone(pytz.timezone("Pacific/Auckland"))    
    
    # approximate arrival time of teleseisms
    surface_wavespeed = 5
    travel_time_sec = great_circ_dist_km / surface_wavespeed
    arrival_time = startLOC + timedelta(seconds=travel_time_sec)
    
    print('='*79)
    print('Event Information')
    print('Origin Time (UTC):',startUTC)
    print('Origin Time (NZ):',startLOC)
    print('Magnitude (Mwc):',event.magnitudes[0].mag)
    print('Latitude:',evnt_lat)
    print('Longitude:',evnt_lon)
    print('Backazimuth:',BAz)
    print('Distance (km):',great_circ_dist_km)
    print('Distance (deg):',approximate_gcd_deg)
    print('Approximate arrival (NZ):',arrival_time)
    print('='*79)
    
    import ipdb;ipdb.set_trace()

        
if __name__ == '__main__':
    event = search_GCMT('2017-298')
    if event:
        fathom_event(event)
    else:
        print('No events found')