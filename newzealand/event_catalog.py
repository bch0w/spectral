"""
The catalog used for initial testing purposes for Pyatoa + Seisflows
26.7.19
"""
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime, Catalog, read_events
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn.header import FDSNNoDataException
from pyatoa.plugins.new_zealand.gather import geonet_mt


def cut_cat(cat, indices=None, method="remove"):
    """
    Obspy's catalog object has no delete, so here it is!

    :type cat: obspy.Catalog
    :param cat: catalog object to cut down
    :type indices: list of ints
    :param indices: indices of the catalog to remove
    :type method: str
    :param method: "remove" to delete 'indices' from 'cat', "keep" to retain
    :rtype: obspy.Catalog
    :return: new catalog object cut from old one
    """
    if indices is None:
        return cat

    new_cat = Catalog()
    if method == "remove":
        events = [e for i, e in enumerate(cat) if i not in indices]
    elif method == "keep":
        events = [e for i, e in enumerate(cat) if i in indices]
    new_cat.extend(events)

    return new_cat


def remove_groupings(cat, sep_km=30.):
    """
    Loop through a catalog and premove spatially grouped events by comparing
    inter-event distances and removing any events that are below some
    separation threshold

    :type cat: obspy.Catalog
    :param cat: catalog to remove grouped events from
    :type sep_km: float
    :param sep_km: desired separation distance between events in km
    :rtype: obspy.Catalog
    :return: new catalog with grouped events removed
    """
    indices_remove = []

    for i, event in enumerate(cat):
        lat = event.preferred_origin().latitude
        lon = event.preferred_origin().longitude

        # Loop through the events that are not this event
        for j, event_check in enumerate(
                [e for j, e in enumerate(cat) if j != i]):
            # Make sure indexing stays the same
            if j >= i:
                j += 1
            # Determine the distance between events
            distance_m, _, _ = gps2dist_azimuth(
                lat1=lat, lon1=lon,
                lat2=event_check.preferred_origin().latitude,
                lon2=event_check.preferred_origin().longitude
            )
            if distance_m < (sep_km * 1E3):
                # Index pairs are the same so avoid counting them twice
                if j in np.unique(indices_remove):
                    continue

                # Pick the event with the larger magnitude
                if (event.preferred_magnitude().mag >
                        event_check.preferred_magnitude().mag):
                    indices_remove.append(j)
                else:
                    indices_remove.append(i)

    print("group removal removing {} events".format(
        len(np.unique(indices_remove))))

    return cut_cat(cat, indices_remove, method="remove")


def spatially_vary(cat, desired_length=None):
    """
    Create a new catalog by picking the next event that is the furthest away 
    from the current event. This way we ensure that our catalog is spatially
    varied and contains all edge events rather than choosing nearest neighbor
    groupings

    NOTE:
    This doesn't really do what I want it to because it still retains grouped
    events after 2 iterations (going back and forth between the two map corners
    for example). What I really want it to do is pick the next event that is the
    furthest away from ALL events in the current catalog, but that would take a 
    lot more checking.

    How I might go about that: 
        -Create a list of distances of the picked event from all events in base
         catalog, sort min to max
        -Create a list of distances of the picked event from all events in the
         new catalog, sort min to max
        -Find the largest shared index between these two lists, corresponding to
         an event with the largest distance from

    :type cat: obspy.Catalog
    :param cat: catalog to spatially vary
    :type desired_length: int
    :param desired_length: desired length of new catalog
    :return:
    """
    if desired_length:
        if len(cat) <= desired_length:
            print("catalog length is less than desired_length")
            return cat
    else:
        desired_length = len(cat)

    # Loop through once to collect latitude/longitude values
    lats, lons = [], []
    for event in cat:
        lats.append(event.preferred_origin().latitude)
        lons.append(event.preferred_origin().longitude)

    # Loop through again to pick values
    varied_indices = []

    # Start at the first index
    i = 0
    varied_indices.append(i)

    # Go until the desired length is hit
    while len(varied_indices) < desired_length:
        # Get the event location
        lat = cat[i].preferred_origin().latitude
        lon = cat[i].preferred_origin().longitude

        # Loop through the events that are not this one
        distances = []
        for j, (lat_chk, lon_chk) in enumerate(
                    [e for j, e in enumerate(zip(lats, lons)) if j != i]):
            # Keep indices consistent
            if j >= i:
                j += 1
            # Keep track of the distances from this event, and their indices
            distance_m, _, _ = gps2dist_azimuth(lat1=lat, lon1=lon,
                                                lat2=lat_chk, lon2=lon_chk
                                                )
            distances.append((j, distance_m))

        # Sort the distances by largest
        largest_distance = sorted(distances, key=lambda t: t[1], reverse=True)

        # Loop through these distances
        for dist in largest_distance:
            # Make sure we haven't already chosen this index
            if dist[0] in varied_indices:
                continue
            print("{} is {:.2f} km from {}".format(
                varied_indices[-1], dist[1] * 1E-3, dist[0]))
            # Append the index and move on
            varied_indices.append(dist[0])
            break

        # Set i to the index of the event we just found
        i = dist[0]

    return varied_indices


def check_moment_tensor(geonet_csv_file, cat):
    """
    Check if the catalog has a moment tensor defined by GeoNet
    Performed based on event ids that are common to both catalogs

    :type geonet_csv_file: str
    :param geonet_csv_file: path to the GeoNet csv file downloaded from
        https://github.com/GeoNet/data
    :type cat: obspy.Catalog
    :param cat: catalog object to check event ids from
    :rtype: obspy.Catalog
    :return: new catalog that only contains events with moment tensors
    """
    # Get all the available event ids in the csv file
    event_ids = []
    with open(geonet_csv_file, 'r') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            event_ids.append(row[0])

    # Loop through the catalog and figure out which event ids are present
    indices_remove = []
    for i, event in enumerate(cat):
        event_id = event.resource_id.id.split('/')[1]
        if event_id not in event_ids:
            indices_remove.append(i)

    print("moment tensor check removed {} events".format(
        len(indices_remove))
    )

    return cut_cat(cat, indices_remove, method="remove")


def append_mt(cat, csv_fid=None):
    """
    Append GeoNet moment tensor information to the events in catalog
    Assumes that all the events have moment tensor information, meaning 
    check_moment_tensor() should be run prior to this

    :type cat: obspy.Catalog
    :param cat: catalog of events to add moment tensor information to
    :type csv_fid: str
    :param csv_fid: csv file containing moment tensor information
    """
    cat_out = Catalog()
    events = []
    for event in cat:
        event_id = event.resource_id.id.split('/')[1]
        event_out, _ = geonet_mt(event_id=event_id, units="dynecm",
                                 event=event, csv_fid=csv_fid)
        events.append(event_out)

    cat_out.extend(events)

    return cat_out


def plot_cat(cat, evnts, tag=""):
    """
    Plot the catalog with a custom title based on the

    :type cat: obspy.Catalog
    :param cat: catalog to plot
    :type evnts: dictionary
    :param evnts: parameter dictionary from function event_trials()
    :type tag: str
    :param tag: tag for saving the plot name
    """
    # Plot the catalog
    evnts['starttime'].precision = 1
    evnts['endtime'].precision = 1

    title = "\n".join([
        "{c} {n} events".format(c=evnts['name'], n=len(cat)),
        "{s} - {e}".format(s=evnts['starttime'], e=evnts['endtime']),
        "M{a} - {b}".format(a=evnts['minmagnitude'], b=evnts['maxmagnitude']),
        "LLC: {latmin},{lonmin} / URC: {latmax},{lonmax}\n".format(
            latmin=evnts['minlatitude'], lonmin=evnts['minlongitude'],
            latmax=evnts['maxlatitude'], lonmax=evnts['maxlongitude'])
    ])

    cat.plot(projection='local', resolution='l', continent_fill_color='w',
             water_fill_color='w', color='date',
             outfile="./{}{}".format(evnts['name'], tag), title=title
             )

    plt.close('all')


def event_trials(choice, csv_file, sep_km=None, desired_length=None):
    """
    For a given set of trial parameters, retrieve an event catalog and then
    slim it down with various check parameters; make sure events are not too 
    closely grouped, and that they contain moment tensor information from GeoNet
    Plot the final catalog, make beachballs for plotting purposes, return
    an .xml catalog

    :type choice: str
    :param choice: choice of parameter trial
    :type csv_file: str
    :param csv_file: path to GeoNet csv file with moment tensor information
    :type desired_length: int
    :param desired_length: desired length of catalog, event removal will try to 
    :return:
    """
    trials = {
        "charlie_trial": {
            "name": "charlie_trial",
            "starttime": UTCDateTime("2010-01-01T00:00:00"),
            "endtime": UTCDateTime("2019-11-01T00:00:00"),
            "minmagnitude": 4.75,
            "maxmagnitude": 6.,
            "minlatitude": -42.5,  # LLC
            "minlongitude": 173.5,  # LLC
            "maxlatitude": -37.25,  # URC
            "maxlongitude": 178.5,  # URC
            "maxdepth": 60.
            },
        "fullscale": {
            "name": "fullscale",
            "starttime": UTCDateTime("2000-01-01T00:00:00"),
            "endtime": UTCDateTime("2019-11-01T00:00:00"),
            "minmagnitude": 4.5,
            "maxmagnitude": 6.,
            "minlatitude": -42.5,  # LLC
            "minlongitude": 173.5,  # LLC
            "maxlatitude": -37.25,  # URC
            "maxlongitude": 178.5,  # URC
            "maxdepth": 60.
          },
        "aspen": {
            "name": "aspen",
            "starttime": UTCDateTime("2003-08-20T00:00:00"),  # GeoNet MT cat
            "endtime": UTCDateTime(),
            "minmagnitude": 4.4,
            "maxmagnitude": 6.,
            "minlatitude": -42.5,  # LLC
            "minlongitude": 173.5,  # LLC
            "maxlatitude": -37.0,  # URC
            "maxlongitude": 178.5,  # URC
            "maxdepth": 60.
          },
        "south": {
            "name": "south",
            "starttime": UTCDateTime("2003-08-20T00:00:00"),  # GeoNet MT cat
            "endtime": UTCDateTime(),
            "minmagnitude": 4.5,
            "maxmagnitude": 6.,
            "minlatitude": -48,  # LLC
            "minlongitude": 165,  # LLC
            "maxlatitude": -40.0,  # URC
            "maxlongitude": 175,  # URC
            "maxdepth": 60.
          }
    }

    evnts = trials[choice]
    
    # Create the catalog and plot the raw event catalog
    original_cat_fid = f"{cat_name}_original.xml"
    if not os.path.exists(original_cat_fid):
        c = Client("GEONET")
        original_cat = c.get_events(starttime=evnts['starttime'],
                                    endtime=evnts['endtime'],
                                    minmagnitude=evnts['minmagnitude'],
                                    maxmagnitude=evnts['maxmagnitude'],
                                    minlatitude=evnts['minlatitude'],
                                    minlongitude=evnts['minlongitude'],
                                    maxlatitude=evnts['maxlatitude'],
                                    maxlongitude=evnts['maxlongitude'],
                                    maxdepth=evnts['maxdepth'],
                                    orderby="magnitude-asc"
                                    )
        print("catalog has {} events".format(len(original_cat)))
        # Delete unncessary attributes to save space
        for attr in ["picks", "amplitudes", "station_magnitudes"]:
            for event in original_cat:
                try:
                    delattr(event, attr)
                except KeyError:
                    continue
        original_cat.write(original_cat_fid, format="QUAKEML")
    else:
        original_cat = read_events(original_cat_fid)

    # Remove if no GeoNet moment tensors
    new_cat = check_moment_tensor(csv_file, original_cat)
    print("catalog has {} events".format(len(new_cat)))

    if sep_km is not None and desired_length is not None:
        cat_out = remove_groupings(new_cat, sep_km=sep_km)
        # Try separation distances until new catalog has desired length
        while True:
            if len(cat_out) == desired_length:
                break
            else:
                print(f"{len(cat_out)} events; " 
                      f"desired: {desired_length}; "
                      f"current separation: {sep_km}")
                sep_km = input("new separation distance?: ")
                if sep_km:
                    sep_km = float(sep_km)

            cat_out = remove_groupings(new_cat, sep_km=sep_km)
            print("catalog has {} events".format(len(cat_out)))
    else:
        cat_out = new_cat

    # add moment tensor information to all events
    cat_out_w_mt = append_mt(cat_out, csv_file)

    cat_out_w_mt.write(f"{cat_name}_w_mt.xml", format="QUAKEML")


def cat_from_event_ids(cat_name, event_ids, csv_file):
    """
    Create an ObsPy catalog with moment tensors from a list of event ids
    """
    c = Client("GEONET")
    cat = Catalog()
    for eventid in event_ids:
        try:
            cat += c.get_events(eventid=eventid)[0]
        except FDSNNoDataException:
            print(f"No event for {eventid}")

    new_cat = check_moment_tensor(csv_file, cat)
    print("catalog has {} events".format(len(new_cat)))

    cat_w_mt = append_mt(new_cat, csv_file)
    cat_w_mt.write(f"{cat_name}_w_mt.xml", format="QUAKEML")



if __name__ == "__main__":
    # Use GeoNet moment tensors
    geonet_path = "geonet/data/moment-tensor/GeoNet_CMT_solutions.csv"
    for path_to in ["/Users/chowbr/Documents/subduction/data",
                    "/Users/Chow/Documents/academic/vuw/data",
                    "/seis/prj/fwi/bchow/data"]:
        if os.path.exists(os.path.join(path_to, geonet_path)):
            csv_file = os.path.join(path_to, geonet_path)
            break
    else:
        csv_file = None
        print("No GeoNet .csv file found")


    # USER PARAMETERS
    cat_name = "south"

    if True:
        # Get the catalog and its name, based on the functions
        event_trials(cat_name, csv_file=csv_file)
    else:
        # Get catalog from list of event ids
        fid = ("/Users/Chow/Documents/academic/vuw/forest/posthoc/"
               "ristau_starred_events.txt")
        event_ids = np.loadtxt(fid, dtype=str)
        cat_from_event_ids(cat_name, event_ids, csv_file)

