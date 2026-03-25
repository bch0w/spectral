"""
TODO Merge this into fcnt2mseed.py

Convert Fairfield node naming convention to chosen station names based on a 
dictionary
"""
from glob import glob
from obspy import read


# For NMSS24 Oct deployment <COMPASS1><BAZ3><ELEMENT1>
dictionary = {
        "100": "S2061",
        "30":  "S2064",
        "101": "S2112",
        "102": "S2113",
        "103": "S3251",
        "104": "S3255",
        "20":  "N2051", 
        "21":  "N2055",
        "22":  "N2012",
        "23":  "N2015",
        "24":  "N2785",
        "25":  "N2786",
        }

for fid in glob("XX.*"):
    net, sta, channel, year, jday = fid.split(".")

    new_sta = dictionary[sta]

    # Change output name
    fid_out = f"{net}.{new_sta}.{channel}..{year}.{jday}"
    print(f"{fid} -> {fid_out}")

    # Change actual stats
    st = read(fid)
    for tr in st:  # Should only be lenght 1
        tr.stats.station = new_sta
    
    st.write(f"name_change/{fid_out}", format="MSEED")


