"""standalone function to generate CMTSOLUTION files for specfem runs
followd the format of Harvard CMTSOLUTION, but replaces centroid moment
entries with more precise GEONET information and MT from John Ristau
"""
import os
import sys
import pandas as pd
sys.path.append('../../modules/')
from getdata import get_moment_tensor, pathnames
from synmod import mt_transform
from obspy.clients.fdsn import Client
from obspy.geodetics import FlinnEngdahl

def get_region(event_list):
    c = Client('GEONET')
    cat = c.get_events(eventid=event_id)
    event = cat[0]
    origin = event.origins[0]
    fe = FlinnEngdahl()
    region = fe.get_region(longitude=origin.longitude,
                  latitude=origin.latitude)

    return event,region

def generate_CMTSOLUTION(event_id):
    """generate CMTSOLUTION file in the format of the Harvard CMT catalog
    -Moment tensor components taken from John Ristaus MT catalog
    -Event information taken from GEONET earthquake catalog
    -CMT information taken from GCMT catalog
    NOTE: template stolen from obspy
    """

    # grab moment tensor information from Ristau's solutions
    MT = get_moment_tensor(event_id=event_id)
    if not MT:
        sys.exit('incorrect event call')
    mt = [MT['Mxx'],MT['Myy'],MT['Mzz'],MT['Mxy'],MT['Mxz'],MT['Myz']]
    mt = [_*(1E20) for _ in mt]
    mt = mt_transform(mt,method='xyz2rtp')
    mrr,mtt,mpp,mrt,mrp,mtp = mt

    event,region = get_region([event_id])
    origin = event.origins[0]
    datetime = origin.time

    # always set zero
    time_shift = 0
    half_duration = 0
    template = (
        "{hypocenter_catalog:>4} {year:4d} {month:02d} {day:02d} {hour:02d} "
        "{minute:02d} {second:05.2f} "
        "{latitude:9.4f} {longitude:9.4f} {depth:5.1f} {mb:.1f} {ms:.1f} "
        "{region}\n"
        "event name:{event_name:>17}\n"
        "time shift:{time_shift:17.4f}\n"
        "half duration:{half_duration:14.4f}\n"
        "latitude:{cmt_latitude:19.4f}\n"
        "longitude:{cmt_longitude:18.4f}\n"
        "depth:{cmt_depth:22.4f}\n"
        "Mrr:{m_rr:24.6E}\n"
        "Mtt:{m_tt:24.6E}\n"
        "Mpp:{m_pp:24.6E}\n"
        "Mrt:{m_rt:24.6E}\n"
        "Mrp:{m_rp:24.6E}\n"
        "Mtp:{m_tp:24.6E}\n"
        )

    template = template.format(hypocenter_catalog="XXXX",
                                year=datetime.year,
                                month=datetime.month,
                                day=datetime.day,
                                hour=datetime.hour,
                                minute=datetime.minute,
                                second=float(datetime.second) +
                                datetime.microsecond / 1E6,
                                latitude=origin.latitude,
                                longitude=origin.longitude,
                                depth=origin.depth / 1000.0,
                                mb=0,
                                ms=0,
                                region=region,
                                event_name=event_id,
                                time_shift=0,
                                half_duration=0,
                                cmt_latitude=origin.latitude,
                                cmt_longitude=origin.longitude,
                                cmt_depth=origin.depth / 1000.0,
                                m_rr=mrr,
                                m_tt=mtt,
                                m_pp=mpp,
                                m_rt=mrt,
                                m_rp=mrp,
                                m_tp=mtp
                                )


    # write to solution file
    filename = (pathnames()['data'] +
                'KUPEDATA/CMTSOLUTIONS/{}CMTSOLUTION'.format(event_id))
    with open(filename,'w') as f:
        f.write(template)
    print(filename)

def mass_get_region(event_list):
    c = Client('GEONET')
    events,regions = [],[]
    for pd_event in event_list:
        event_id = pd_event['event_id']
        cat = c.get_events(eventid=event_id)
        event = cat[0]
        origin = event.origins[0]
        fe = FlinnEngdahl()
        region = fe.get_region(longitude=origin.longitude,
                      latitude=origin.latitude)
        
        events.append(event)
        regions.append(region)
        
    return events,regions

def generate_CMTSOLUTION_from_tomCat(event_id=None):
    """generate CMTSOLUTION file from tomCat
    28.3 haven't tested it yet
    """
    tomCat_path = pathnames()['data'] + 'tomCat/tomCat'
    tomCat = pd.read_pickle(tomCat_path)
    template = (
        "{hypocenter_catalog:>4} {year:4d} {month:02d} {day:02d} {hour:02d} "
        "{minute:02d} {second:05.2f} "
        "{latitude:9.4f} {longitude:9.4f} {depth:5.1f} {mb:.1f} {ms:.1f} "
        "{region}\n"
        "event name:{event_name:>17}\n"
        "time shift:{time_shift:17.4f}\n"
        "half duration:{half_duration:14.4f}\n"
        "latitude:{cmt_latitude:19.4f}\n"
        "longitude:{cmt_longitude:18.4f}\n"
        "depth:{cmt_depth:22.4f}\n"
        "Mrr:{m_rr:24.6E}\n"
        "Mtt:{m_tt:24.6E}\n"
        "Mpp:{m_pp:24.6E}\n"
        "Mrt:{m_rt:24.6E}\n"
        "Mrp:{m_rp:24.6E}\n"
        "Mtp:{m_tp:24.6E}\n"
        )
    
    event_list = []
    # single event
    if event_id:
        event = tomCat.loc[tomCat['event_id'] == event_id].iloc[0]
        if event.empty:
            print("{} empty".format(event_id))
            return
        else:
            event_list.append(event)
    # mass process all events in tomCat
    else:
        for index,event in tomCat.iterrows():
            event_list.append(event)
    print(len(event_list),'events to convert')
    
    # get regions beforehand to save time
    _,regions = mass_get_region(event_list)
    
    for event,region in zip(event_list,regions):
        # parse tomCat
        event_id = event['event_id']
        datetime = event['datetime']
        latitude = event['latitude']
        longitude = event['longitude']
        depth = event['depth']
        mrr = event['m_rr']
        mtt = event['m_tt']
        mpp = event['m_pp']
        mrt = event['m_rt']
        mrp = event['m_rp']
        mtp = event['m_tp']
        settemplate = template.format(hypocenter_catalog="XXXX",
                                        year=datetime.year,
                                        month=datetime.month,
                                        day=datetime.day,
                                        hour=datetime.hour,
                                        minute=datetime.minute,
                                        second=float(datetime.second) +
                                        datetime.microsecond / 1E6,
                                        latitude=latitude,
                                        longitude=longitude,
                                        depth=depth,
                                        mb=0,
                                        ms=0,
                                        region=region,
                                        event_name=event_id,
                                        time_shift=0,
                                        half_duration=0,
                                        cmt_latitude=latitude,
                                        cmt_longitude=longitude,
                                        cmt_depth=depth,
                                        m_rr=mrr,
                                        m_tt=mtt,
                                        m_pp=mpp,
                                        m_rt=mrt,
                                        m_rp=mrp,
                                        m_tp=mtp
                                        )

        # write to solution file
        filename = (pathnames()['data'] +
                    'KUPEDATA/CMTSOLUTIONS/{}CMTSOLUTION'.format(event_id))
        print(os.path.basename(filename))
        with open(filename,'w') as f:
            f.write(settemplate)
    
if __name__ == "__main__":
    # generate_CMTSOLUTION(sys.argv[1])
    generate_CMTSOLUTION_from_tomCat()
