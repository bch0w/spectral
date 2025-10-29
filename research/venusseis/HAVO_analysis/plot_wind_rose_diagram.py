"""
Thanks, ChatGPT!
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

from obspy import read

# -------------------------------------------------------
# Example synthetic data (replace these with your own)
# -------------------------------------------------------
fid = sys.argv[1]
st = read(fid).merge()
wind_dir_deg = st[0].data
n = len(wind_dir_deg)

start_time = st[0].stats.starttime.datetime
wind_time = np.array([start_time + timedelta(hours=i) for i in range(n)])

# -------------------------------------------------------
# Define function to categorize time of day
# -------------------------------------------------------
def time_of_day(hour):
    """times given in UTC, shift to HST which is UTC-10"""
    # hour = hour - 10 % 24
    if 0 <= hour < 6:
        return 0  # Morning
    elif 6 <= hour < 12:
        return 1  # Afternoon
    elif 12 <= hour < 18:
        return 2  # Evening
    else:
        return 3  # Night

# Map each timestamp to a time-of-day category (0–3)
hours = np.array([t.hour for t in wind_time])
period_indices = np.array([time_of_day(h) for h in hours])

period_labels = ["UTC 00-06", "UTC 06-12", "UTC 12-18", "UTC 18-24"]
num_periods = len(period_labels)

# -------------------------------------------------------
# Define direction bins
# -------------------------------------------------------
num_bins = 64
bin_edges = np.linspace(0, 360, num_bins + 1)
bin_centers = np.deg2rad((bin_edges[:-1] + bin_edges[1:]) / 2)

# -------------------------------------------------------
# Count occurrences per direction per time period
# -------------------------------------------------------
counts = np.zeros((num_periods, num_bins))

for p in range(num_periods):
    mask = (period_indices == p)
    hist, _ = np.histogram(wind_dir_deg[mask], bins=bin_edges)
    counts[p, :] = hist

# Normalize counts (optional)
counts = counts / counts.sum()

# -------------------------------------------------------
# Plot: stacked wind rose colored by time of day
# -------------------------------------------------------
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 8))

bottom = np.zeros(num_bins)
colors = plt.cm.plasma_r(np.linspace(0, 1, num_periods))

for p in range(num_periods):
    bars = ax.bar(bin_centers, counts[p], width=(2 * np.pi / num_bins),
                  bottom=bottom, color=colors[p],
                  edgecolor='black', linewidth=0.5, alpha=0.8,
                  label=period_labels[p])
    bottom += counts[p]

# -------------------------------------------------------
# Style the plot
# -------------------------------------------------------
ax.set_theta_zero_location("N")   # North on top
ax.set_theta_direction(-1)        # Clockwise
ax.set_title(os.path.basename(fid), va='bottom')
ax.set_yticklabels([])            # Hide radial labels
ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))

plt.tight_layout()
plt.show()




