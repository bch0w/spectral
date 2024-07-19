"""
Declustering source stations, need to get a count of how many stations each
source has data for to maximize the amount of data included in the inversion
"""
import os
import numpy as np
from glob import glob


# Filenames
stapairs = "/Users/chow/Repos/spectral/research/alaska/AKATOM/station_pairs.txt"
stations = "/Users/chow/Repos/spectral/research/alaska/AKATOM/nalaska/DATA/STATIONS"
sources = "/Users/chow/Repos/spectral/research/alaska/AKATOM/nalaska/DATA/FORCESOLUTIONS_ALL"

# Get source names
source_names = glob(f"{sources}/*")
source_names = [os.path.basename(_)[14:] for _ in source_names]
source_names = [_.replace("_", ".") for _ in source_names]

# Get station names
station_names = []
with open(stations, "r") as f:
    for line in f.readlines():
        sta, net, *_ = line.split()
        station_names.append(f"{net}.{sta}")

# Create station pairs dictionary
station_pairs = {}
src = None
record = False
with open(stapairs, "r") as f:
    lines = f.readlines()
for line in lines:
    # Denotes a source station
    if line.startswith("RAY") or line.startswith("LOV"):
        wave, phase, src = line.strip().split()

        # We only want available 'HYP' data
        if phase == "ELL":
            record = False
        elif src not in source_names:
            record = False
        else:
            record = True
            print(line)
            station_pairs[src] = []

    # This is a source station we want to record data for
    elif record:
        station = line.strip()
        if station in station_names:
            print(line)
            station_pairs[src].append(station)



