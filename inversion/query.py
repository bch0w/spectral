"""
Providing an event id will return the classification of the event by querying
a CSV file that contains that information
"""
import sys
import csv


def query_single(event_id):
    """
    Return the stats of a single event
    """
    with open("./events.csv", "r") as f:
        events = csv.reader(f, delimiter=",")
        for e in events:
            if e[0] == event_id:
                print(f"\nId:   {e[0]}\n"
                      f"Time: {e[1]}\n"
                      f"Mag:  {e[2]}\n"
                      f"Dpt:  {e[3]}\n"
                      f"Lat:  {e[4]}\n"
                      f"Lon:  {e[5]}\n"
                      f"Qty:  {e[6]}\n")
                break
        else:
            print(f"{event_id} not found")


def query_list(event_list):
    """
    User-provided list returns quality for each event
    """
    quality = [-1] * len(event_list)
    with open("./events.csv", "r") as f:
        events = csv.reader(f, delimiter=",")
        for e in events:
            if e[0] in event_list:
                idx = event_list.index(e[0])
                quality[idx] = int(e[-1])
    
    # sort the two lists together
    quality, event_list = zip(*sorted(zip(quality, event_list)))

    for q, e in zip(quality, event_list):
        print(f"{e:<11}: {q:>3}")


if __name__ == "__main__":
    # query_single(sys.argv[1])
    query_list()

