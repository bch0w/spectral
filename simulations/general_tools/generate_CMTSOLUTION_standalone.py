"""standalone function to generate CMTSOLUTION files for specfem runs
following the format of Harvard CMTSOLUTION, but replaces centroid moment
entries with more precise GEONET information and MT from John Ristau
"""
import sys
import csv
from obspy.clients.fdsn import Client
from obspy.geodetics import FlinnEngdahl


def get_moment_tensor(event_id, csv_file):
    """
    gets moment tensor as array from geonet CSV file
    """
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            if i == 0:
                tags = row
            if event_id == row[0]:
                values = []
                for t, v in zip(tags, row):
                    if (t == "Date") or (t == "PublicID"):
                        values.append(v)
                    else:
                        values.append(float(v))

                MT = dict(zip(tags, values))
                return MT


def mt_transform(mt, method):
    """
    transform moment tensor between xyz and rtp coordinates
    acceptable mt formats:
        [m11,m22,m33,m12,m13,m23]
        [mxx,myy,mzz,mxy,mxz,myz]
        [mrr,mtt,mpp,mrt,mrp,mtp]
    :type mt: dict
    :param mt: moment tensor in format above
    :type method: str
    :param method: type of conversion, "rtp2xyz" or "xyz2rtp"
    """
    if method == "xyz2rtp":
        if "m_xx" not in mt.keys():
            print("for xyz2rtp, dict must have keys in xyz")
        m_rr = mt["m_zz"]
        m_tt = mt["m_xx"]
        m_pp = mt["m_yy"]
        m_rt = mt["m_xz"]
        m_rp = -1 * mt["m_yz"]
        m_tp = -1 * mt["m_xy"]
        return {"m_rr": m_rr, "m_tt": m_tt, "m_pp": m_pp, "m_rt": m_rt,
                "m_rp": m_rp, "m_tp": m_tp}

    if method == "rtp2xyz":
        if "m_tt" not in mt.keys():
            print("for rtp2xyz, dict must have keys in rtp")
        m_xx = mt["m_tt"]
        m_yy = mt["m_pp"]
        m_zz = mt["m_rr"]
        m_xy = -1 * mt["m_tp"]
        m_xz = mt["m_rt"]
        m_yz = -1 * mt["m_rp"]
        return {"m_xx": m_xx, "m_yy": m_yy, "m_zz": m_zz, "m_xy": m_xy,
                "m_xz": m_xz, "m_yz": m_yz}
    else:
        print("Invalid transformation method, xyz2rtp or rtp2xyz")
        return None


def get_event_and_region(event_or_id):
    """
    get region for a more complete looking CMTSOLUTION file
    :param event_or_id:
    :return:
    """
    if isinstance(event_or_id, str):
        c = Client('GEONET')
        cat = c.get_events(eventid=event_or_id)
        event = cat[0]
    else:
        event = event_or_id
    origin = event.origins[0]
    fe = FlinnEngdahl()
    region = fe.get_region(longitude=origin.longitude, latitude=origin.latitude)

    return event, region


def generate_CMTSOLUTION(event_or_id, csv_file, output_file):
    """generate CMTSOLUTION file in the format of the Harvard CMT catalog
    -Moment tensor components taken from John Ristaus MT catalog
    -Event information taken from GEONET earthquake catalog
    -CMT information taken from GCMT catalog
    NOTE: template stolen from obspy
    """
    if not isinstance(event_or_id, str):
        event_id = event_or_id.resource_id.id.split('/')[1]
    else:
        event_id = event_or_id

    # grab moment tensor information from Ristau's solutions
    MT = get_moment_tensor(event_id=event_id, csv_file=csv_file)
    if not MT:
        print("{} not found".format(event_id))
        return
    mt = {"m_xx": MT['Mxx'], "m_yy": MT['Myy'], "m_zz": MT['Mzz'],
          "m_xy": MT['Mxy'], "m_xz": MT['Mxz'], "m_yz": MT['Myz']
          }
    for key in mt:
        mt[key] *= 1E20
    mt = mt_transform(mt, method='xyz2rtp')

    event, region = get_event_and_region(event_or_id)
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
    
    template = template.format(hypocenter_catalog="XXXX", year=datetime.year,
                               month=datetime.month, day=datetime.day,
                               hour=datetime.hour, minute=datetime.minute,
                               second=float(datetime.second) +
                                                datetime.microsecond / 1E6,
                               latitude=origin.latitude,
                               longitude=origin.longitude,
                               depth=origin.depth / 1000.0, mb=0, ms=0,
                               region=region, event_name=event_id,
                               time_shift=time_shift,
                               half_duration=half_duration,
                               cmt_latitude=origin.latitude,
                               cmt_longitude=origin.longitude,
                               cmt_depth=origin.depth / 1000.0, m_rr=mt["m_rr"],
                               m_tt=mt["m_tt"], m_pp=mt["m_pp"],
                               m_rt=mt["m_rt"], m_rp=mt["m_rp"], m_tp=mt["m_tp"]
                               )

    # write to solution file
    filename = output_file
    with open(filename, 'w') as f:
        f.write(template)
    print(filename)


if __name__ == "__main__":
    # ======================== set parameters ==================================
    eventid = "2016p355601"
    # csv_file = "path/to/GeoNet_CMT_solutions.csv"
    csv_file = "/Users/chowbr/Documents/subduction/data/GEONET/data/moment-tensor/GeoNet_CMT_solutions.csv"
    output_path = "/Users/chowbr/Documents/subduction/data/KUPEDATA/CMTSOLUTIONS/{}CMTSOLUTION".format(eventid)
    # ======================== set parameters ==================================

    generate_CMTSOLUTION(eventid, csv_file, output_path)
