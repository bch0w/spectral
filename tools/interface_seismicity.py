"""
An earthquake catalog, a plate interface model, and user-defined parameters walk
into a bar...
Sort through the earthquake catalog based on vertical distance to the plate
interface. Brute force because I'm too lazy to do this nicely.
Earthquake file should be in the .csv output from GeoNet QuakeSearch
"""
import sys
import numpy as np


# REQUIRED USER PARAMETER
drange = 1000.  # allowable distance to interface


def myround(x, base=5, choice='near'):
    """
    Round value x to nearest base, round 'up','down' or to 'near'est base

    :type x: float
    :param x: value to be rounded
    :type base: int
    :param base: nearest integer to be rounded to
    :type choice: str
    :param choice: method of rounding, 'up', 'down' or 'near'
    :rtype roundout: int
    :return: rounded value
    """
    if choice == 'near':
        roundout = int(base * round(float(x)/base))
    elif choice == 'down':
        roundout = int(base * np.floor(float(x)/base))
    elif choice == 'up':
        roundout = int(base * np.ceil(float(x)/base))

    return roundout


# Open the input data
assert(".csv" in sys.argv[1]), "Input file must be a .csv"
eqs = np.genfromtxt(sys.argv[1], delimiter=",", dtype=str)
intfc = np.load("williams_hikurangi_interface.npy")

# A bit wonky but we're going to convert all coordinates to ints to truncate
# floating points at 3 dec. places to avoid floating point rounding errors
intfc = (intfc * 1E3).astype(int)

# Calculate grid spacing on the interface file
xvals = np.unique(intfc[:, 0])
yvals = np.unique(intfc[:, 1])
dx = xvals[1] - xvals[0]
dy = yvals[1] - yvals[0]

# Grab relevant information from earthquake csv file
eqs = eqs[1:, (2, 4, 5, 6, 7)]
print(f"{len(eqs)} earthquakes in catalog")
with open(f"{sys.argv[1]}_filtered.txt", "w") as f:
    for i, eq in enumerate(eqs):
        time = eq[0]
        lon, lat, mag, zkm = eq[1:].astype(float)

        # Do the same integer conversion to truncate floating points
        lon_ = int(myround(x=lon*1E3, base=dx))
        lat_ = int(myround(x=lat*1E3, base=dx))
        zkm_ = int(zkm * 1E3)

        # Find corresponding index and depth value
        idx = np.where((intfc[:, 0] == lat_) & (intfc[:, 1] == lon_))[0]
        if idx.size == 1:
            # Finally, check if the earthquake depth falls near the interface
            # Interface values are negative but eq depths are positive
            plate = intfc[idx][0][2]
            depth_range = (abs(plate) - drange,  abs(plate) + drange)
            if (zkm_ > depth_range[0]) and (zkm_ < depth_range[1]):
                print(f"{i:0>3}/{len(eqs)}")
                f.write(f"{lon}\t{lat}\n")
