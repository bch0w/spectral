"""
Earthquake getter script based on station or array locations. Allow the user
to define epicentral rings around a geographic location that are linked to
magnitude which empirically accounts for amplitude attenuation via geometric
spreading and intrinsic attenuation.

For example, within a ring of 20km, collect earthquakes that have any
magnitude. But from 20--200km, cap the minimum magnitude at M3. And from
200-1000km, cap the minimum magnitude at M5 etc. This allows gathering a
catalog of events that does not intrinsically throw out low-magnitude
earthquakes, and that attempts to maximize the number of events with good SNR
recorded by a station or array over a given time period.

Multiple locations (e.g., each station of an array) may be given; rings are
queried around each location and the resulting catalogs are merged with
duplicate events removed.

.. rubric::

    # Edit the PARAMETERS section below, then run
    python acquire_staggered_quakeml.py
"""
from obspy import UTCDateTime, Catalog
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy.geodetics import kilometer2degrees


# =============================== PARAMETERS ==================================
# FDSN client with an event service (EarthScope no longer serves events)
CLIENT = "USGS"
CHOICE = "AK2"

# Distances define a radius in km within which we allow the corresponding
# magnitude or larger (i.e., the magnitude defines the lower limit). Rings are
# annular, i.e., the 100 km entry covers 20--100 km. The last entry (global)
# defines the lowest magnitude acceptable from the largest radius out to the
# antipode. A magnitude of None means any magnitude is accepted
MAG_DIST_DICT = {
    50: 2,
    500: 3,
    5000: 5,
    "global": 7,
}
PLOT = True

if CHOICE == "HI":
    # Time bounds of the catalog
    STARTTIME = UTCDateTime("2025-03-25T00:00:00")
    ENDTIME = UTCDateTime("2025-05-04T00:00:00")

    # Optional length of time chunks [days] for each query, to avoid hitting
    # server-side event limits for long durations. Set None to query all at once
    CHUNK_DAYS = 365

    # List of (latitude, longitude) reference points, e.g., a single station, an
    # array centroid, or every station in an array
    LOCATIONS = [
        (19.335887, -155.301148),  # MITD
    ]

    # Optional depth bounds [km], set None to ignore
    MINDEPTH = None
    MAXDEPTH = None

    # Output QuakeML file
    OUTPUT = "ishtar-hi_staggered_catalog.xml"

elif CHOICE == "AK1":
    # Time bounds of the catalog
    STARTTIME = UTCDateTime("2025-07-09T00:00:00")
    ENDTIME = UTCDateTime("2025-07-23T00:00:00")

    # Optional length of time chunks [days] for each query, to avoid hitting
    # server-side event limits for long durations. Set None to query all at once
    CHUNK_DAYS = 365

    # List of (latitude, longitude) reference points, e.g., a single station, an
    # array centroid, or every station in an array
    LOCATIONS = [
        (65.82507, -149.54317),  # H23K
        # (64.9602, -148.2319),  # MDM
    ]

    # Optional depth bounds [km], set None to ignore
    MINDEPTH = None
    MAXDEPTH = None

    # Output QuakeML file
    OUTPUT = "ishtar-ak1_staggered_catalog.xml"

elif CHOICE == "AK2":
    # Time bounds of the catalog
    STARTTIME = UTCDateTime("2025-09-02T00:00:00")
    ENDTIME = UTCDateTime("2025-09-29T00:00:00")

    # Optional length of time chunks [days] for each query, to avoid hitting
    # server-side event limits for long durations. Set None to query all at once
    CHUNK_DAYS = 365

    # List of (latitude, longitude) reference points, e.g., a single station, an
    # array centroid, or every station in an array
    LOCATIONS = [
        (62.96992, -145.46989),  # PAX
    ]

    # Optional depth bounds [km], set None to ignore
    MINDEPTH = None
    MAXDEPTH = None

    # Output QuakeML file
    OUTPUT = "ishtar-ak2_staggered_catalog.xml"

# =============================================================================


def build_rings(mag_dist_dict):
    """
    Convert the distance-magnitude dictionary into a list of annular rings
    sorted by increasing radius.

    :type mag_dist_dict: dict
    :param mag_dist_dict: keys are outer ring radii in km (or 'global' for
        the outermost ring extending to the antipode), values are minimum
        magnitudes (or None for no magnitude limit)
    :rtype: list of tuple
    :return: list of (minradius_deg, maxradius_deg, minmagnitude) tuples
    """
    radii = sorted(k for k in mag_dist_dict if k != "global")
    rings = []
    inner = 0.
    for radius in radii:
        outer = kilometer2degrees(radius)
        rings.append((inner, outer, mag_dist_dict[radius]))
        inner = outer
    if "global" in mag_dist_dict:
        rings.append((inner, 180., mag_dist_dict["global"]))

    return rings


def time_chunks(starttime, endtime, chunk_days=None):
    """
    Split a time range into consecutive chunks.

    :type starttime: obspy.UTCDateTime
    :param starttime: start of the time range
    :type endtime: obspy.UTCDateTime
    :param endtime: end of the time range
    :type chunk_days: float or None
    :param chunk_days: length of each chunk in days. If None, the whole time
        range is returned as a single chunk
    :rtype: list of tuple
    :return: list of (chunk_start, chunk_end) tuples
    """
    if not chunk_days:
        return [(starttime, endtime)]

    chunks = []
    t = starttime
    while t < endtime:
        t_end = min(t + chunk_days * 86400, endtime)
        chunks.append((t, t_end))
        t = t_end

    return chunks


def query_ring(client, lat, lon, minradius, maxradius, minmagnitude,
               starttime, endtime, chunk_days=None, mindepth=None,
               maxdepth=None):
    """
    Query events within a single annular ring around a location.

    :type client: obspy.clients.fdsn.Client
    :param client: FDSN client used to query events
    :type lat: float
    :param lat: latitude of the ring center
    :type lon: float
    :param lon: longitude of the ring center
    :type minradius: float
    :param minradius: inner ring radius in degrees
    :type maxradius: float
    :param maxradius: outer ring radius in degrees
    :type minmagnitude: float or None
    :param minmagnitude: minimum allowable magnitude, None for no limit
    :type starttime: obspy.UTCDateTime
    :param starttime: start time of the query
    :type endtime: obspy.UTCDateTime
    :param endtime: end time of the query
    :type chunk_days: float or None
    :param chunk_days: optional query chunk length in days
    :type mindepth: float or None
    :param mindepth: optional minimum event depth in km
    :type maxdepth: float or None
    :param maxdepth: optional maximum event depth in km
    :rtype: obspy.Catalog
    :return: events found within the ring
    """
    cat = Catalog()
    kwargs = dict(latitude=lat, longitude=lon, minradius=minradius,
                  maxradius=maxradius, minmagnitude=minmagnitude,
                  mindepth=mindepth, maxdepth=maxdepth)
    # Drop unset parameters so that the server applies no constraint
    kwargs = {k: v for k, v in kwargs.items() if v is not None}

    for t0, t1 in time_chunks(starttime, endtime, chunk_days):
        try:
            cat += client.get_events(starttime=t0, endtime=t1, **kwargs)
        except FDSNNoDataException:
            continue

    return cat


def remove_duplicates(cat):
    """
    Remove duplicate events from a catalog based on their resource IDs.

    :type cat: obspy.Catalog
    :param cat: catalog that may contain duplicate events
    :rtype: obspy.Catalog
    :return: catalog with unique events, sorted by origin time
    """
    unique = {}
    for event in cat:
        unique.setdefault(str(event.resource_id), event)

    events = sorted(
        unique.values(),
        key=lambda e: (e.preferred_origin() or e.origins[0]).time
    )
    return Catalog(events=events)


def acquire_staggered_catalog(client, locations, mag_dist_dict, starttime,
                              endtime, chunk_days=None, mindepth=None,
                              maxdepth=None):
    """
    Gather a catalog of events using distance-dependent minimum magnitudes
    around one or more reference locations.

    :type client: obspy.clients.fdsn.Client
    :param client: FDSN client used to query events
    :type locations: list of tuple
    :param locations: (latitude, longitude) reference points
    :type mag_dist_dict: dict
    :param mag_dist_dict: ring radii in km (and 'global') mapped to minimum
        magnitudes, see `build_rings`
    :type starttime: obspy.UTCDateTime
    :param starttime: start time of the catalog
    :type endtime: obspy.UTCDateTime
    :param endtime: end time of the catalog
    :type chunk_days: float or None
    :param chunk_days: optional query chunk length in days
    :type mindepth: float or None
    :param mindepth: optional minimum event depth in km
    :type maxdepth: float or None
    :param maxdepth: optional maximum event depth in km
    :rtype: obspy.Catalog
    :return: merged catalog of unique events
    """
    rings = build_rings(mag_dist_dict)
    cat = Catalog()
    for lat, lon in locations:
        print(f"location: {lat:.4f}, {lon:.4f}")
        for minradius, maxradius, minmag in rings:
            ring_cat = query_ring(client, lat, lon, minradius, maxradius,
                                  minmag, starttime, endtime, chunk_days,
                                  mindepth, maxdepth)
            mag_str = "any" if minmag is None else f">= M{minmag}"
            print(f"\t{minradius:7.3f}--{maxradius:7.3f} deg, {mag_str}: "
                  f"{len(ring_cat)} events")
            cat += ring_cat

    return remove_duplicates(cat)


if __name__ == "__main__":
    client = Client(CLIENT)
    cat = acquire_staggered_catalog(
        client=client, locations=LOCATIONS, mag_dist_dict=MAG_DIST_DICT,
        starttime=STARTTIME, endtime=ENDTIME, chunk_days=CHUNK_DAYS,
        mindepth=MINDEPTH, maxdepth=MAXDEPTH
    )
    print(f"{len(cat)} unique events, writing to '{OUTPUT}'")
    cat.write(OUTPUT, format="QUAKEML")
    if PLOT:
        cat.plot(proj="ortho")
