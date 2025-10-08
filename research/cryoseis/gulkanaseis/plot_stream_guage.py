import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.imaging.util import ObsPyAutoDateFormatter
from dateutil.rrule import MINUTELY, SECONDLY
from matplotlib.dates import (date2num, AutoDateLocator, AutoDateFormatter, 
                              DateFormatter, num2date)


def _set_xaxis_obspy_dates(ax, ticklabels_small=True):
    """
    Set Formatter/Locator of x-Axis to use ObsPyAutoDateFormatter and do some
    other tweaking.

    In contrast to normal matplotlib ``AutoDateFormatter`` e.g. shows full
    timestamp on first tick when zoomed in so far that matplotlib would only
    show hours or minutes on all ticks (making it impossible to tell the date
    from the axis labels) and also shows full timestamp in matplotlib figures
    info line (mouse-over info of current cursor position).

    :type ax: :class:`matplotlib.axes.Axes`
    :rtype: None
    """
    ax.xaxis_date()
    locator = AutoDateLocator(minticks=3, maxticks=50)
    locator.intervald[MINUTELY] = [1, 2, 5, 10, 15, 30]
    locator.intervald[SECONDLY] = [1, 2, 5, 10, 15, 30]
    ax.xaxis.set_major_formatter(ObsPyAutoDateFormatter(locator))
    ax.xaxis.set_major_locator(locator)
    if ticklabels_small:
        plt.setp(ax.get_xticklabels(), fontsize='small')

# SET 
path = "/Users/chow/Work/research/gulkanaseis24/data/USGS_data/phelan_creek_stream_guage_2024-09-07_to_2024-09-14.txt"
data = np.loadtxt(path, skiprows=28, usecols=[2,4], delimiter="\t", dtype=str)
times, height_ft = data.T  # time in AK local

# Convert arrays
times = np.array([UTCDateTime(_).datetime for _ in times])
height_m = np.array([_ * 0.3048 for _ in height_ft.astype(float)])

plt.plot(times, height_m, lw=1)
_set_xaxis_obspy_dates(plt.gca())

# Tmarks

plt.title("USGS Phelan Cr. Stream Guage (63.24027778, -145.4680556)")
plt.ylabel("Height [m]")
plt.xlabel("Time [AKDT]")
plt.grid()

plt.show()

