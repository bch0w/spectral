"""
Standalone function to generate CMTSOLUTION files from an ObsPY catalog
Only usable for the New Zealand problem, hardcoding for GeoNet moment tensors
Expects that the events have been gathered using the 'event_catalog.py' script
"""
import os
import sys
from obspy import read_events
from obspy.geodetics import FlinnEngdahl


def get_region(event):
    """
    Get region for a more complete looking CMTSOLUTION file

    :type event: obspy.event
    :param event: event
    :rtype: str
    :return: Flinn Engdahl region based on lat lon
    """
    origin = event.origins[0]
    fe = FlinnEngdahl()
    region = fe.get_region(longitude=origin.longitude, latitude=origin.latitude)

    return region


def generate_cmtsolutions(catalog, convert=1E7, path="./"):
    """
    Generate CMTSOLUTION file in the format of the Harvard CMT catalog
    -Moment tensor components taken from John Ristaus MT catalog
    -Event information taken from GEONET earthquake catalog
    -CMT information taken from GCMT catalog

    NOTE: template stolen from obspy

    :type catalog: obspy.core.event.catalog.Catalog
    :param catalog: Catalog that contains focal mechanisms to use for CMTs
    :type convert: float
    :param convert: to convert the units of the moment tensors. Use-case:
        ObsPy moment tensors take moment tensor components in units of N*m, but
        CMTSOLUTIONS expect units in dyne*cm, which means conversion by 1E7 is
        required to get the correct format. Defaults to 1.
    :type path: str
    :param path: path to save the CMTSOLUTION files to
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
        assert(hasattr(event, "focal_mechanisms")), "Event {event_id} has no " \
                                                    "focal mechanism"

        # Moment tensors have already been converted to GCMT format
        mt = event.preferred_focal_mechanism().moment_tensor.tensor

        # Get some auxiliary information
        region = get_region(event)
        origin = event.origins[0]
        datetime = origin.time

        # Get magnitude information for quick reference
        if hasattr(event, "magnitudes"):
            mag = event.preferred_magnitude().mag
        else:
            mag = 0

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
                                   depth=origin.depth / 1000.0, mb=0, ms=mag,
                                   region=region, event_name=event_id,
                                   time_shift=time_shift,
                                   half_duration=half_duration,
                                   cmt_latitude=origin.latitude,
                                   cmt_longitude=origin.longitude,
                                   cmt_depth=origin.depth / 1000.0,
                                   m_rr=mt.m_rr * convert, 
                                   m_tt=mt.m_tt * convert, 
                                   m_pp=mt.m_pp * convert,
                                   m_rt=mt.m_rt * convert, 
                                   m_rp=mt.m_rp * convert, 
                                   m_tp=mt.m_tp * convert
                                   )

        # Write to solution file
        filename = os.path.join(path, "CMTSOLUTION_{}".format(event_id))
        with open(filename, 'w') as f:
            f.write(data_out)

        print(filename)


if __name__ == "__main__":
    cat = read_events(sys.argv[1])
    generate_cmtsolutions(catalog=cat, convert=1, path="./")
