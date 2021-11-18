"""
Prep a catalog to run through Carl's decluster matlab scripts. It requires files
to be in the format:

event_id,time,amg,depth_km,lat,lon
"""
import numpy as np
from obspy import read_events, Catalog

lat_min = -42.5
lat_max = -37.0
lon_min = 173.0
lon_max = 178.5
depth_min_km = 60
mag_min = 4.5
mag_max = 6.

cat = read_events("/Users/Chow/Documents/academic/vuw/forest/utils/prep/posthoc/nz_north_cat_253.xml")
forest = np.loadtxt("/Users/Chow/Documents/academic/vuw/forest/utils/prep/posthoc/prepfiles/decluster_60_ids.txt", dtype=str)
alpha30 = np.loadtxt("/Users/Chow/Documents/academic/vuw/forest/utils/prep/posthoc/alpha_final_30event/final_posthoc_event_ids.txt", dtype=str)
skip = np.append(forest, alpha30)
print(f"SKIP: {len(skip)}")

cat_out = []
for event in cat:
    eid = event.resource_id.id.split("/")[1]
    print(eid)
    lat = event.preferred_origin().latitude
    lon = event.preferred_origin().longitude
    dpt = event.preferred_origin().depth * 1E-3
    mag = event.preferred_magnitude().mag

    if (eid in skip) or (lat <= lat_min) or (lat >= lat_max) or \
            (lon <= lon_min) or (lon >= lon_max) or (dpt > depth_min_km) or \
            (mag < mag_min) or (mag >= mag_max):
        print(f"{eid} skipped for something")
        continue
    else:
        cat_out.append(event)

cat_write = Catalog()
cat_write.extend(cat_out)
cat_write.write(f"new_catalog_{len(cat_out)}.xml", format="quakeml")

