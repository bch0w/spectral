"""
Standalone function to generate CMTSOLUTION files from an obspy catalog
Only usable for the New Zealand problem, hardcoding for GeoNet moment tensors
"""
import os
import sys
import csv
from obspy.clients.fdsn import Client
from obspy.geodetics import FlinnEngdahl
from obspy.imaging.beachball import beachball


def get_moment_tensor(event_id, csv_file):
    """
    gets moment tensor as array from geonet CSV file
    :type event_id: str
    :param event_id: GeoNet event id, e.g. 2018p130600
    :type csv_file: str
    :param csv_file: path to the GeoNet CSV file containing moment tensors
        https://github.com/GeoNet/data/blob/master/moment-tensor/
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

                mt = dict(zip(tags, values))
                return mt


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


def get_event_and_region(event):
    """
    get region for a more complete looking CMTSOLUTION file
    :type event: obspy.event
    :param event: event
    :return:
    """
    origin = event.origins[0]
    fe = FlinnEngdahl()
    region = fe.get_region(longitude=origin.longitude, latitude=origin.latitude)

    return event, region


def make_beachball(mt, outfile="./beachball.png"):
    """
    This is a convenient time to plot beachballs for use in plotting
    """
    width = 2.6E4
    facecolor = "r"

    beach_input = [mt['m_rr'], mt['m_tt'], mt['m_pp'], 
                   mt['m_rt'], mt['m_rp'], mt['m_tp']]

    b = beachball(beach_input, facecolor="r", outfile=outfile)
    

def generate_cmtsolutions(catalog, csv_file, make_bball=False, pathout="./"):
    """
    Generate CMTSOLUTION file in the format of the Harvard CMT catalog
    -Moment tensor components taken from John Ristaus MT catalog
    -Event information taken from GEONET earthquake catalog
    -CMT information taken from GCMT catalog
    NOTE: template stolen from obspy
    """
    template = (
        "{hypocenter_cat:>4} {year:4d} {month:02d} {day:02d} {hour:02d} "
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

    for event in catalog:
        event_id = event.resource_id.id.split('/')[1]

        # grab moment tensor information from Ristau's solutions
        moten = get_moment_tensor(event_id=event_id, csv_file=csv_file)
        if not moten:
            print("{} not found".format(event_id))
            continue

        # Transform Ristau's moment tensor to follow Harvarad CMT catalog format
        mt = {"m_xx": moten['Mxx'], "m_yy": moten['Myy'], "m_zz": moten['Mzz'],
              "m_xy": moten['Mxy'], "m_xz": moten['Mxz'], "m_yz": moten['Myz']
              }
        # Convert to dyne*cm
        for key in mt:
            mt[key] *= 1E20
        mt = mt_transform(mt, method='xyz2rtp')
        
        # make beachballs
        if make_bball:
            make_beachball(mt, outfile=f"{event_id}_beachball.png")

        # Get some auxiliary information
        event, region = get_event_and_region(event)
        origin = event.origins[0]
        datetime = origin.time

        # Always set these zero for simulations
        time_shift = 0
        half_duration = 0

        # Apply information to the template
        data_out = template.format(hypocenter_cat="XXXX",
                                   year=datetime.year,
                                   month=datetime.month, day=datetime.day,
                                   hour=datetime.hour, minute=datetime.minute,
                                   second=(float(datetime.second) +
                                           datetime.microsecond / 1E6),
                                   latitude=origin.latitude,
                                   longitude=origin.longitude,
                                   depth=origin.depth / 1000.0, mb=0, ms=0,
                                   region=region, event_name=event_id,
                                   time_shift=time_shift,
                                   half_duration=half_duration,
                                   cmt_latitude=origin.latitude,
                                   cmt_longitude=origin.longitude,
                                   cmt_depth=origin.depth / 1000.0,
                                   m_rr=mt["m_rr"], m_tt=mt["m_tt"],
                                   m_pp=mt["m_pp"], m_rt=mt["m_rt"],
                                   m_rp=mt["m_rp"], m_tp=mt["m_tp"]
                                   )

        # Write to solution file
        filename = os.path.join(pathout, "CMTSOLUTION_{}".format(event_id))
        with open(filename, 'w') as f:
            f.write(data_out)
        print(filename)


if __name__ == "__main__":
    pass
