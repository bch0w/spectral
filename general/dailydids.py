"""
A script for simple journal entries which tags time, location and tasks
"""
import os
import sys
import datetime

DATETIME_STR = "%Y-%m-%d %H:%M:%S"


def print_last(fid, n=1):
    """
    Print the last `n` entries

    :type fid: str
    :param fid: File containing the journal entries
    :type n: int
    :param n: number of entries to print, defaults to 1
    """
    with open(fid, "r") as f:
        lines = f.readlines()
        entries = 0
        for i, line in enumerate(reversed(lines)):
            if "IN  :" in line or "OUT :" in line:
                entries += 1
                if entries >= n:
                    break
        print("\n")
        print("".join(lines[-1*i-1:]))


def check_last_entry(fid):
    """
    See if the last entry is IN or OUT, if IN create new entry
    """
    with open(fid, "r") as f:
        lines = f.readlines()
        # read file from the end, backwards
        for line in reversed(lines):
            if "IN  :" in line:
                # Get the last time in entry    
                time_in = ":".join(line.split(":")[1:]).strip()
                time_in = datetime.datetime.strptime(time_in, DATETIME_STR)
                break
            elif "OUT :" in line:
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
    Create a journal entry for what to do today
    """
    NOW = datetime.datetime.now()

    in_or_out = input("'in' or 'out'?: ").upper()
    if in_or_out not in ["IN", "OUT"]:
        return

    if in_or_out == "IN":
        check_last_entry(fid=fid)

    with open(fid, "a") as f:
        # Write the time in or out
        entry_time = input("time in/out? (leave blank if now, or give minutes "
                           "before now): ")
        if entry_time:
            entry_time = NOW - datetime.timedelta(minutes=int(entry_time))
        else:
            entry_time = NOW

        f.write(f"{in_or_out:<4}: {entry_time:{DATETIME_STR}}\n")

        # Determine where you're working
        if in_or_out == "IN":
            where = input("where? 0-UAF, 1-HOME, 3-OTHER: ")
            where_list = ["UAF", "HOME"]
            try:
                where = where_list[int(where)]
            except IndexError:
                where = input("other?: ").upper()
            f.write(f"WORK: {where}\n")
            # Write down the to do list for today
            while True:
                to_do = input("to do today?: ")
                if not to_do:
                    check = input("is that all? ([y]/n): ")
                    if check != "n":
                        break
                f.write(f"TODO: {to_do} \n")
        # Leaving for the day, only track time and tasks accomplshed
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
    # Ensure the file exists
    # TERN UAF
    txtfile = "/home/bchow/work/repos/spectral/general/dailydids.txt"
    if not os.path.exists(txtfile):
        raise FileNotFoundError(f"No such file {txtfile}")

    # Either check the last entry or write a new entry
    try:
        sys.argv[1] = "last"
        try:
            n = int(sys.argv[2])
        except IndexError:
            n = 1
        print_last(fid=txtfile, n=n)
    except IndexError:
        daily_did(fid=txtfile)



