from obspy import UTCDateTime

with open("dailydids.txt", "r") as f:
    lines = f.readlines()


day1 = 287

days = []
try:
    for line in lines:
        if "IN  :" in line:
            time_in = UTCDateTime(line.split(":")[1])
        elif "OUT :" in line:
            time_out = UTCDateTime(line.split(":")[1])
            elapsed = (time_out - time_in) / 3600
            if time_out.year == 2018:
                day = time_out.julday - day1
            elif time_out.year > 2018:
                # Rest of 2018 + Each additional year + Remainder of year
                day = (365 - 287) + 365 * (time_out.year - 2019) + time_out.julday
            days.append((day, elapsed))
except Exception as e:
    print(e)
    import ipdb;ipdb.set_trace()


import ipdb;ipdb.set_trace()
