"""
Script to fetch GeoNet earthquake information, write out relevant information
into a version controlled CSV file for parsing, and also provide the ability
to output event information into CMTSOLUTION files to be fed into specfem, as
well as into .XML files of manageable size for quick event fetching
"""
from obspy import UTCDateTime


def fetch_events(starttime="2017-01-01T00:00:00", output_fid=None):
    """
    fetch events using the obspy FDSN interface
    :param starttime: catalog start time
    :param output_fid: file path for catalog write
    :return cat: obspy catalog object
    """
    from obspy.clients.fdsn import Client

    c = Client("GEONET")
    cat = c.get_events(starttime=UTCDateTime(starttime), endtime=UTCDateTime(),
                       minmagnitude=4.5, minlatitude=-42.5007,
                       maxlatitude=-36.9488, minlongitude=172.9998,
                       maxlongitude=179.5077
                       )

    if output_fid:
        cat.write(output_fid, format="QUAKEML")

    return cat


def write_catalog_to_csv(cat, output_fid="./tomcat.csv"):
    """
    as the name says
    :param cat:
    :return:
    """
    import csv
    import ipdb;ipdb.set_trace()

    with open(output_fid, 'w') as csvfile:
        fieldnames = ["Event ID", "Origin Time", "Magnitude", "Latitude",
                      "Longitude", "Depth"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for event in cat:
            event_id = event.resource_id.id.split('/')[1]
            time = str(event.preferred_origin().time)
            magnitude = "{}{:.2f}".format(
                event.preferred_magnitude().magnitude_type,
                event.preferred_magnitude().mag)
            latitude = event.preferred_origin().latitude
            longitude = event.preferred_origin().longitude
            depth = "{:.2f}".format(event.preferred_origin().depth*1E-3)
            writer.writerow({"Event ID": event_id, "Origin Time": time,
                             "Magnitude": magnitude, "Latitude": latitude,
                             "Longitude": longitude, "Depth": depth}
                            )


def output_as_cmtsolutions(cat):
    """
    write a catalog as CMTSOLUTIONS, check if files already exist
    :param cat:
    :return:
    """
    import sys
    sys.path.append('../../simulations/general_tools')
    from os.path import exists
    from importlib import reload
    import generate_CMTSOLUTION_standalone
    reload(generate_CMTSOLUTION_standalone)
    import generate_CMTSOLUTION_standalone as genny

    csv_file = ("/Users/chowbr/Documents/subduction/data/GEONET/data/"
                "moment-tensor/GeoNet_CMT_solutions.csv")
    output_path = ("/Users/chowbr/Documents/subduction/data/KUPEDATA/" 
                   "CMTSOLUTIONS/{}CMTSOLUTION")

    for event in cat:
        event_id = event.resource_id.id.split('/')[1]
        event_output = output_path.format(event_id)
        if exists(event_output):
            continue
        genny.generate_CMTSOLUTION(event_or_id=event, csv_file=csv_file,
                                   output_file=event_output)


if __name__ == "__main__":
    today = str(UTCDateTime().date)
    cat = fetch_events(output_fid="{}_EVENT_CATALOG.xml".format(today))
    write_catalog_to_csv(cat, output_fid="{}_EVENT_CATALOG.csv".format(today))
    output_as_cmtsolutions(cat)
