import os
import datetime

datetime_str = "%Y-%m-%d %I:%M:%S"

def check_last_entry(fid):
    """
    see if the last entry is IN or OUT, if IN create new entry
    """
    with open(fid, "r") as f:
        lines = f.readlines()
        # read file from the end, backwards
        for line in reversed(lines):
            if "IN" in line:
                # Get the last time in entry    
                time_in = ":".join(line.split(":")[1:]).strip()
                time_in = datetime.datetime.strptime(time_in, datetime_str)
                break
            elif "OUT" in line:
                return

    # Will only be accessed if last entry was IN
    print(f"No matching OUT entry for last IN entry...")
    with open(fid, "a") as f:
        # Query the time left
        time_delta = input(f"time out in minutes after {time_in}?: ")
        time_out = time_in + datetime.timedelta(minutes=int(time_delta))
        f.write(f"{'OUT':<4}: {time_out}\n")
        while True:
            done = input("accomplished that day?: ")
            if not done:
                check = input("is that all? ([y]/n): ")
                if check != "n":
                    break
            f.write(f"DONE: {done} \n")
        f.write("\n")
    print("\n")


def daily_did(fid):
    """
    create a journal entry for what to do today
    """
    now = datetime.datetime.now()
    in_or_out = input("'in' or 'out'?: ")
    if in_or_out not in ["in", "out"]:
        return
    in_or_out = in_or_out.upper()
    if in_or_out == "IN":
        check_last_entry(fid=fid)

    with open(fid, "a") as f:
        # Write the time in or out
        entry_time = input("time in/out? (leave blank if now, or give minutes before now): ")
        if entry_time:
            entry_time = now - datetime.timedelta(minutes=int(entry_time))
        else:
            entry_time = now

        f.write(f"{in_or_out:<4}: {entry_time:{datetime_str}}\n")

        # Determine where you're working
        if in_or_out == "IN":
            where = input("where? 0-VUW, 1-GNS, 2-HOME, 3-OTHER: ")
            where_list = ["VUW", "GNS", "HOME"]
            try:
                where = where_list[int(where)]
            except IndexError:
                where = input("other?: ").upper()
            f.write(f"WORK: {where}\n")

        # coming in for the day
        if in_or_out == "IN":
            while True:
                to_do = input("to do today?: ")
                if not to_do:
                    check = input("is that all? ([y]/n): ")
                    if check != "n":
                        break
                f.write(f"TODO: {to_do} \n")
        # leaving for the day
        if in_or_out == "OUT":
            while True:
                done = input("accomplished today?: ")
                if not done:
                    check = input("is that all? ([y]/n): ")
                    if check != "n":
                        break
                f.write(f"DONE: {done} \n")
        f.write("\n")
               
if __name__ == "__main__":
    txtfile = "/Users/Chow/Documents/academic/vuw/packages/spectral/general/dailydids.txt"
    if not os.path.exists(txtfile):
        raise FileNotFoundError(f"No such file {txtfile}")
    daily_did(fid=txtfile)



