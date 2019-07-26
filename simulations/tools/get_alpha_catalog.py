"""
The catalog used for initial testing purposes for Pyatoa + Seisflows
26.7.19
"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import catalog_to_cmtsolutions

catalog_name = "alpha_trial"

c = Client("GEONET")
starttime = UTCDateTime("2012-01-01T00:00:00")  # GeoNet changed how they label
endtime = UTCDateTime("2019-01-01T00:00")       # event ID's in 2012
minmagnitude = 5.  # Large enough events to get good observation signals
maxmagnitude = 6.  # Small enough events to avoid finite-fault effects
minlatitude = -41.5  # LLC to ignore Kaikoura aftershock sequence
minlongitude = 173.  # LLC
maxlatitude = -37.3  # URC low enough to exclude 2016 Te Araroa sequence,
maxlongitude = 179.  # which was poorly located due to its location

# This catalog excludes Kaikoura and Te Araroa Sequences
cat = c.get_events(starttime=starttime, endtime=endtime,
                   minmagnitude=minmagnitude, maxmagnitude=maxmagnitude,
                   minlatitude=minlatitude, minlongitude=minlongitude,
                   maxlatitude=maxlatitude, maxlongitude=maxlongitude
                   )

# Plot the catalog
starttime.precision = 1
endtime.precision = 1
title = "\n".join(["{c} {n} events".format(c=catalog_name, n=len(cat)),
                   "{s} - {e}".format(s=starttime, e=endtime),
                   "M{a} - {b}".format(a=minmagnitude, b=maxmagnitude),
                   "LLC: {latmin},{lonmin} / URC: {latmax},{lonmax}\n".format(
                       latmin=minlatitude, lonmin=minlongitude,
                       latmax=maxlatitude, lonmax=maxlongitude)]
                  )
cat.plot(projection='local', resolution='l', continent_fill_color='w',
         water_fill_color='w', color='date',
         outfile="./{}".format(catalog_name), title=title
         )

# Write out the Quakeml file
cat.write("{}.xml".format(catalog_name), format="QUAKEML")

# Create all the CMTSOLUTIONS
csv_file = "/Users/chowbr/Documents/subduction/data/GEONET/data/" \
           "moment-tensor/GeoNet_CMT_solutions.csv"

# Four events will not have GeoNet moment tensors, leaving 26 trial events
generate_cmtsolutions(cat, csv_file)


