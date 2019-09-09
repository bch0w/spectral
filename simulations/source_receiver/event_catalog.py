"""
The catalog used for initial testing purposes for Pyatoa + Seisflows
26.7.19
"""
import csv
import numpy as np
from obspy import UTCDateTime, Catalog
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth

import sys
sys.path.append(".")
import catalog_to_cmtsolutions


def cut_cat(cat, indices=[], method="remove"):
    """
    Obspy's catalog object has no delete, so here it is!
    :return:
    """
    new_cat = Catalog()
    if method == "remove":
        events = [e for i, e in enumerate(cat) if i not in indices]
    elif method == "keep":
        events = [e for i, e in enumerate(cat) if i in indices]
    new_cat.extend(events)

    return new_cat


def mesh_domain(cat, lat_min, lon_min, lat_max, lon_max, distance_km=50.):
    """
    Keep events away from the mesh boundaries to avoid spurious reflections

    :param cat:
    :param lat_min:
    :param lon_min:
    :param lat_max:
    :return:
    """
    raise NotImplementedError


def remove_groupings(cat, sep_km=30.):
    """
    Loop through a catalog and premove spatially grouped events
    :param cat:
    :return:
    """
    indices_remove = []

    # Stop when the number of events to remove hits a desired length
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
    Loop through a catalog and pick the event that is the furthest away from
    the current event. This way we ensure that our catalog is spatially varied.

    NOTE:
    This doesn't really do what I want it to because it still retains grouped
    events after 2 iterations (going back and forth between the two map corners
    for example). What I really want it to do is pick the next event that is the
    furthest away from ALL events in the catalog, but that would take a lot of
    checking.

    :param cat:
    :param desired_length:
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
    :return:
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


def plot_cat(cat, evnts, tag=""):
    """
    Plot the catalog with a custom title based on the
    :param cat:
    :param evnts:
    :param tag:
    :return:
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
             outfile="./{}_{}".format(evnts['name'], tag), title=title
             )


def alpha_trial():
    """
    Alpha trial gives 26 events. Used for initial checkerboard testing in
    late August, early September, 2019.
    2012 starttime to avoid the event ID name change from GeoNet
    Top corner picked to avoid Te Araroa sequence
    Bottom corner chosen to avoid numerous Kaikoura aftershocks
    :return:
    """
    evnts = {
        "name": "alpha_trial",
        "starttime": UTCDateTime("2012-01-01T00:00:00"),
        "endtime": UTCDateTime("2019-01-01T00:00"),
        "minmagnitude": 5.,
        "maxmagnitude": 6.,
        "minlatitude": -41.5,  # LLC
        "minlongitude": 173.,  # LLC
        "maxlatitude": -37.3,  # URC
        "maxlongitude": 179.,  # URC
         }

    # Create the catalog and plot the raw event catalog
    c = Client("GEONET")
    cat = c.get_events(starttime=evnts['starttime'], endtime=evnts['endtime'],
                       minmagnitude=evnts['minmagnitude'],
                       maxmagnitude=evnts['maxmagnitude'],
                       minlatitude=evnts['minlatitude'],
                       minlongitude=evnts['minlongitude'],
                       maxlatitude=evnts['maxlatitude'],
                       maxlongitude=evnts['maxlongitude']
                       )

    plot_cat(cat, evnts, tag="_raw_{}".format(len(cat)))

    return cat, evnts["name"]


def beta_trial(csv_file, desired_length):
    """
    An attempt to get more events than the alpha_trial, and to spatially
    distribute them by removing grouped events. 30km separation at group removal
    stage gives 30 events
    :return:
    """
    evnts = {
        "name": "beta_trial",
        "starttime": UTCDateTime("2012-01-01T00:00:00"),
        "endtime": UTCDateTime("2019-01-01T00:00"),
        "minmagnitude": 4.75,
        "maxmagnitude": 6.,
        "minlatitude": -42.5,  # LLC
        "minlongitude": 173.5,  # LLC
        "maxlatitude": -37.25,  # URC
        "maxlongitude": 178.5,  # URC
         }

    # Create the catalog and plot the raw event catalog
    c = Client("GEONET")
    original_cat = c.get_events(starttime=evnts['starttime'],
                                endtime=evnts['endtime'],
                                minmagnitude=evnts['minmagnitude'],
                                maxmagnitude=evnts['maxmagnitude'],
                                minlatitude=evnts['minlatitude'],
                                minlongitude=evnts['minlongitude'],
                                maxlatitude=evnts['maxlatitude'],
                                maxlongitude=evnts['maxlongitude']
                                )
    print("catalog has {} events".format(len(original_cat)))

    # Remove if no GeoNet moment tensors
    new_cat = check_moment_tensor(csv_file, original_cat)
    print("catalog has {} events".format(len(new_cat)))

    # Remove grouped events
    sep_km = 30.
    new_cat = remove_groupings(new_cat, sep_km=sep_km)
    print("catalog has {} events".format(len(new_cat)))

    plot_cat(new_cat, evnts)

    return new_cat, evnts["name"]


if __name__ == "__main__":
    # Use GeoNet moment tensors
    csv_file = ("/Users/chowbr/Documents/subduction/data/GEONET/data/" 
                "moment-tensor/GeoNet_CMT_solutions.csv")

    desired_catalog_length = 30

    # Get the catalog and its name, based on the functions
    cat, cat_name = beta_trial(csv_file, desired_catalog_length)

    # Write to an XML file
    cat.write("{cat_name}.xml", format="QUAKEML")

    # Write the catalog to CMTSOLUTION files required by Specfem3D,
    catalog_to_cmtsolutions.generate_cmtsolutions(cat, csv_file)


