
"""
Parse dailydids and make a plot
"""
import matplotlib.pyplot as plt

from obspy import UTCDateTime


points = []
in_check = 0
lines = open("./dailydids.txt").readlines()

for i, line in enumerate(lines):
    if line[:3] == "IN ":
        time_in = UTCDateTime(line[6:].strip())
        in_check += 1
    elif line[:3] == "OUT":
        time_out = UTCDateTime(line[6:].strip())
        points.append((time_in, time_out))
        in_check = 0
    elif in_check > 1:
        print("multiple in statements with no out statement")
        a=1/0

import ipdb;ipdb.set_trace()

        
