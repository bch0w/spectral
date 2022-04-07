"""
OhRahRah

A graphical illustration of factors that might affect the visibility of the
northern lights (aurora borealis). 
"""
from time import time
import pandas as pd
from urllib.request import urlopen
import matplotlib.pyplot as plt


def get_data_url(url):
    """
    :type url: str
    :param url: url to get data from 
    """
    file = urlopen(url)
    lines = []
    for line in file:
        lines.append(line.decode("utf-8"))

    print(f"retrieved {len(lines)} lines of text from url: {url}")

    return lines

def noaa_mag_7day():
    """
    Get NOAA SWPC data (7day, 1 minute sampling rate) for solar wind 
    (Bx, By, Bz)
    """
    url = "https://services.swpc.noaa.gov/products/solar-wind/mag-7-day.json"
    lines = get_data_url(url)

    from IPython import embed; embed()


def noaa_plasma_7day():
    """
    Get NOAA SWPC data (7day, 1 minute sampling rate) for solar wind 
    (Bx, By, Bz)
    """
    url = "https://services.swpc.noaa.gov/products/solar-wind/plasma-7-day.json"
    lines = get_data_url(url)

    from IPython import embed


def gfz_kp_data():
    """
    GFZ Potsdam produces geomagnetic planetary three-hour index Kp and 
    associated geomagnetic indices as well as relevant solar indices. Pull this
    data from the web and return for plotting 

    :type lines: list
    :param lines: list of lines read from GFZ text file
    :rtype vals: pandas.DataFrame
    :return vals: dataframe of values listed in the text file
    """
    url = ("https://www-app3.gfz-potsdam.de/kp_index/"
           "Kp_ap_Ap_SN_F107_nowcast.txt")
    lines = get_data_url(url)

    vals = {}
    # Grab the header info
    for line in lines:
        if "YYY" in line:
            # Strip the leading '#'
            header = line[1:].strip().split()
            vals = {key: [] for key in header}
        # Grab formatting information
        elif "iii" in line:
            formatter = line[1:].strip().split()

    # Append data to output dictionary
    for line in lines:
        # Skip comments
        if line[0] == "#":
            continue
        else:
            line = line.strip().split()
            for key, val, fmt in zip(header, line, formatter):
                # Want dates as strings so we can concatenate them
                if key in ["YYY", "MM", "DD"]:
                    fmt_ = str
                elif "i" in fmt:
                    fmt_ = int
                elif "f" in fmt:
                    fmt_ = float
                vals[key].append(fmt_(val))

    # Work with Pandas to simplify dealing with this much data
    df = pd.DataFrame(vals)

    # Generate start and endtimes for the dataframe
    time_col = []
    rows, cols = df.shape
    starts = df["YYY"] + "-" + df["MM"] + "-" + df["DD"] + " 00:00:00Z"
    ends = df["YYY"] + "-" + df["MM"] + "-" + df["DD"] + " 11:59:59Z"
    df["starttimes"] = pd.to_datetime(starts)
    df["endtimes"] = pd.to_datetime(ends)

    return df



if __name__ == "__main__":
    # df_kp = gfz_kp_data()
    df_mag = noaa_mag_7day()
    df_pls = noaa_plasma_7day()




