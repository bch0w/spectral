"""
NOAA emails me Kp index forecasts but they're in UTC and hard to interpret 
in plain text. Just plot them as a hist in a new time zone for quick ingestion.
"""
from time import time
import matplotlib.pyplot as plt

text = """NOAA Kp index forecast 17 Feb - 19 Feb
             Feb 17    Feb 18    Feb 19
00-03UT        2         2         2
03-06UT        1         2         2
06-09UT        1         1         3
09-12UT        1         1         3
12-15UT        1         1         3
15-18UT        1         1         4
18-21UT        2         2         4
21-00UT        2         1         4"""

lines = text.split("\n")
title = lines[0]
dates = lines[1].strip().split("   ")
times_utc, day_a, day_b, day_c = [], [], [], []
for line in lines[2:]:
    t, a, b, c = line.strip().split()
    times_utc.append(t)
    day_a.append(int(a))
    day_b.append(int(b))
    day_c.append(int(c))

assert(times_utc[0] == "00-03UT")  # quick check

days = day_a + day_b + day_c
start_time = 15  # convert UTC 00 to local time zone, military time
time = list(range(start_time, start_time + 24, 3))
for i, t in enumerate(time[:]):
    if t >= 24:
        time[i] -= 24
# Setting up the cyclical time axis
times = time * 3
times_label = [f"{t:0>2}" for t in times]
xaxis = range(len(times))

plt.plot(xaxis, days, "r.-")

# Break up days
i = 0
for x, time in enumerate(times):
    if time == 0:
        plt.axvline(x, c="k")
        plt.text(s=dates[i], x=x, y=9)
        i += 1

plt.xticks(xaxis, times_label, rotation="vertical")
plt.xlabel("Time AKST [hour]")
plt.ylabel("NOAA predicted Kp Index")
plt.title(title)
plt.ylim([0,10])
plt.grid(alpha=0.2)

save_fid = title.replace(" ", "_")
save_fid = save_fid.replace("-", "to")
plt.savefig(save_fid)
plt.show()
