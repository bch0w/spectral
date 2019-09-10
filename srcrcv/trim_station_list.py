"""
Remove groupings of stations in temporary network data
"""
import numpy as np
from obspy.geodetics import gps2dist_azimuth


manual_remove = [
    # TVZ LOWER GROUP
    ('NZ', 'COVZ'),
    ('NZ', 'WPVZ'),
    ('NZ', 'FWVZ'),
    ('NZ', 'MAVZ'),
    ('NZ', 'WHVZ'),
    # TVZ UPPER GROUP
    ('NZ', 'NTVZ'),
    ('NZ', 'TMVZ'),
    ('NZ', 'OTVZ'),
]
# Easier to keep these then remove all
manual_keep = []
    # BANNISTER X1
    ('X1', 'HD02'),
    ('X1', 'HD29'),
('X1', 'HD18'),
('X1', 'HD55'),



]

station_file = "./MASTER_STATION_LIST"
stations = np.genfromtxt(station_file, dtype='str')
networks = ['X1', 'X2', 'NZ']

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

# Visual


