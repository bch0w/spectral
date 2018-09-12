"""standalone function to generate CMTSOLUTION files for specfem runs
following the format of Harvard CMTSOLUTION, but replaces centroid moment
entries with more precise GEONET information and MT from John Ristau
"""
import sys
import csv
from obspy.clients.fdsn import Client
from obspy.geodetics import FlinnEngdahl

def get_moment_tensor(event_id,csv_file):
    """gets moment tensor as array from geonet CSV file"""
    with open(csv_file,'r') as f:
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

def mt_transform(mt,method):
    """transform moment tensor between xyz and rtp coordinates
    acceptable mt formats:
        [m11,m22,m33,m12,m13,m23]
        [mxx,myy,mzz,mxy,mxz,myz]
        [mrr,mtt,mpp,mrt,mrp,mtp]
    :type mt: list
    :param mt: moment tensor in format above
    :type method: str
    :param method: type of conversion, "rtp2xyz" or "xyz2rtp"
    """
    if method == 'xyz2rtp':
        m_xx,m_yy,m_zz,m_xy,m_xz,m_yz = mt
        m_rr = m_zz
        m_tt = m_xx
        m_pp = m_yy
        m_rt = m_xz
        m_rp = -m_yz
        m_tp = -m_xy
        return [m_rr,m_tt,m_pp,m_rt,m_rp,m_tp]
    if method == 'rtp2xyz':
        m_rr,m_tt,m_pp,m_rt,m_rp,m_tp = mt
        m_xx = m_tt
        m_yy = m_pp
        m_zz = m_rr
        m_xy = -m_tp
        m_xz = m_rt
        m_yz = -m_rp
        return [m_xx,m_yy,m_zz,m_xy,m_xz,m_yz]
    else:
        print("Invalid transformation method")

def get_event_and_region(event_id):
    c = Client('GEONET')
    cat = c.get_events(eventid=event_id)
    event = cat[0]
    origin = event.origins[0]
    fe = FlinnEngdahl()
    region = fe.get_region(longitude=origin.longitude,
                  latitude=origin.latitude)

    return event,region

def generate_CMTSOLUTION(event_id,csv_file,output_file):
    """generate CMTSOLUTION file in the format of the Harvard CMT catalog
    -Moment tensor components taken from John Ristaus MT catalog
    -Event information taken from GEONET earthquake catalog
    -CMT information taken from GCMT catalog
    NOTE: template stolen from obspy
    """

    # grab moment tensor information from Ristau's solutions
    MT = get_moment_tensor(event_id=event_id,csv_file=csv_file)
    if not MT:
        sys.exit('incorrect event call')
    mt = [MT['Mxx'],MT['Myy'],MT['Mzz'],MT['Mxy'],MT['Mxz'],MT['Myz']]
    mt = [_*(1E20) for _ in mt]
    mt = mt_transform(mt,method='xyz2rtp')
    mrr,mtt,mpp,mrt,mrp,mtp = mt

    event,region = get_event_and_region(event_id)
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
    filename = output_file
    with open(filename,'w') as f:
        f.write(template)
    print(filename)

if __name__ == "__main__":
    # ======================== set parameters ==================================
    eventid = "2014p240655"
    csv_file = "path/to/GeoNet_CMT_solutions.csv"
    output_path = "path/to/OUTPUT_CMTSOLUTION"
    # ======================== set parameters ==================================
    generate_CMTSOLUTION(eventid,csv_file,output_path)
