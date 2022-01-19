from obspy import UTCDateTime

daily_dids = "./dailydids.txt"
lines = open(daily_dids, "r").readlines()

date_convert = {"jan":1, "feb":2, "mar":3, "apr":4, "may":5, "jun":6, "jul":7,
                "aug":8, "sep":9, "oct":10, "nov":11, "dec":12
                }

dates, starts, ends = [],[],[]
for line in lines:
    try:
        if "IN" in line:
            header = line.split()
            if len(header) == 6:
                day, month, day_date, at, loc, hours = header
                start = hours.split("/")[0][1:]
                end = hours.split("/")[1][1:]
                date = UTCDateTime("2010-{month}-{day}".format(
                    month=date_convert[month.lower()], day=day_date)
                )
                dates.append(date)
                starts.append(start)
                ends.append(end)
            else:
                print(header, len(header))
    except Exception as e:
        import ipdb;ipdb.set_trace()
        continue
