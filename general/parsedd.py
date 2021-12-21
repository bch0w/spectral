
"""
Parse dailydids and make a plot
"""
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

points = []
in_check = 0
lines = open("./vuw_dailydids.txt").readlines()
dt_fmt = "%Y-%m-%d %H:%M:%S.%f"
dt_fmt_2 = "%Y-%m-%d %H:%M:%S"  # some dates change fmt and drop the decimal

for i, line in enumerate(lines):
    if line[:3] == "IN ":
        try:
            time_in = dt.datetime.strptime(line[6:].strip(), dt_fmt)
        except ValueError:
            time_in = dt.datetime.strptime(line[6:].strip(), dt_fmt_2)
        in_check += 1
    elif line[:3] == "OUT":
        try:
            time_out = dt.datetime.strptime(line[6:].strip(), dt_fmt)
        except ValueError:
            time_out = dt.datetime.strptime(line[6:].strip(), dt_fmt_2)
        points.append((time_in, time_out))
        in_check = 0
    elif in_check > 1:
        # print("multiple in statements with no out statement, skip")
        continue


# Want to plot hours on the y-axis, dates on the x-axis
ax = plt.subplot()
twax = plt.twinx()
ax.yaxis.set_major_locator(mdates.HourLocator(byhour=range(1)))
# ax.yaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))

time_worked, date_in = [], []
for point in points:
    time_in, time_out = point
    time_in_time = time_in.time().hour + time_in.time().minute / 60.
    time_out_time = time_out.time().hour + time_out.time().minute / 60.
    # Check if times make sense, if not then assume non-military time error
    if time_in_time > time_out_time:
        time_out_time += 12
    delta = time_out_time - time_in_time
    if delta < 4:
        c = "r"
    elif 4 <= delta < 6:
        c = "y"
    elif 6 <= delta <= 8:
        c = "orange"
    else:
        c = "g"
    ax.plot([time_in.date(), time_in.date()], 
            [time_in_time, time_out_time], c=c, alpha=.8, zorder=7)
    time_worked.append(abs(time_out_time - time_in_time))
    date_in.append(time_in.date())

# twax.scatter(date_in, time_worked, c="k", marker="o", s=.5, alpha=0.5)
time_worked_avg = moving_average(time_worked, n=7)
twax.plot(date_in[3:-3], time_worked_avg, "kx-", lw=1, alpha=0.75, markersize=.5, 
          zorder=5)

ax.invert_yaxis()
ax.set_ylim([23, 6])
plt.gcf().autofmt_xdate()

# important dates
imp_dates = {"AGU19": "2019-12-13",
             "AGU20": "2020-12-01",
             "NZ LOCKDOWN": "2020-03-25",
             "HAND IN": "2021-04-30",
             "DEFENSE": "2021-08-12"
             }
for key, value in imp_dates.items():
    date_ = dt.datetime.strptime(value, "%Y-%m-%d")
    plt.axvline(date_, 0, 1, c="k", lw=1)
    ax.text(date_, 7, key, fontsize=8)


plt.xlabel("Date")
ax.set_ylabel("Time of day")
twax.set_ylabel("time worked (1-week rolling average) [h]")
plt.show()
