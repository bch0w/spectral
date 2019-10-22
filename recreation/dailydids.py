import datetime

def daily_did(fid):
    """
    create a journal entry for what to do today
    """
    now = datetime.datetime.now()
    in_or_out = input("'in' or 'out'?: ")
    if in_or_out not in ["in", "out"]:
        return
    in_or_out = in_or_out.upper()

    with open(fid, "a") as f:
        # Write the time in or out
        entry_time = input("time in/out? (leave blank if now, or give minutes before now): ")
        if entry_time:
            entry_time = now - datetime.timedelta(minutes=int(entry_time))
        else:
            entry_time = now
        f.write(f"{in_or_out:<4}: {entry_time}\n")

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
    txtfile = "/Users/chowbr/Documents/subduction/spectral/recreation/dailydids.txt"
    daily_did(fid=txtfile)



