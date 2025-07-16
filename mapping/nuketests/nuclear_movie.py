"""
Generate a series of still frames showing nuclear tests on a map by date
"""
import pandas as pd
import numpy as np
from calendar import monthrange


# Read data and parse out key info
df = pd.read_csv("/Users/chow/Repos/spectral/mapping/nuketests/"
                 "nuclear_explosions.csv")

# Format date to get year, month, day
dates = df.date_long
years = np.array([int(f"{_}"[:4]) for _ in dates])
months = np.array([int(f"{_}"[4:6]) for _ in dates])
days = np.array([int(f"{_}"[7:]) for _ in dates])

# Get locations
lats = df.latitude.to_numpy()
lons = df.longitude.to_numpy()
depths = df.depth.to_numpy()

# Assign colors to each of the countries
country = df.country.to_numpy()
colordict = {"CHINA": "yellow", 
             "FRANCE": "deepskyblue",
             "INDIA": "orange", 
             "PAKIST": "forestgreen",
             "UK": "white",
             "USA": "blue", 
             "USSR": "red"
             }
colors = np.array([colordict[_] for _ in country])

# Get information for size
mb = df.magnitude_body.to_numpy()
yield_ = df.yield_upper.to_numpy()

# Loop throgh all available dates
for year in range(years.min(), years.max() + 1, 1):
    for month in range(1, 12 + 1, 1):
        # Figure out how many days in this particular month 
        start, end = monthrange(year, month)
        for day in range(start, end + 1, 1):
            long_date = int(f"{year}{month:0>2}{day:0>2}")
            idx = df.index[df.date_long == long_date].to_list()
            if idx:


