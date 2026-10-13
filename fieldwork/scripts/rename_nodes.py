"""
Standalone renaming script to change node naming conventions
"""
import sys
from glob import glob
from obspy import read

# USER PARAMETERS
network_dict = {
        "XX": "XX"
        }
station_dict = {
        "101": "HC01",
        "102": "HC02",
        "103": "HC03",
        "104": "HC04",
        "105": "HC05",
        "106": "HC06",
        "107": "HC07",
        "108": "HC08",
        "109": "HC09",
        "110": "HC10",
        "111": "HC11",
        "112": "HC12",
        "113": "HC13",
        "115": "HC15",
        "116": "HC16",
        "201": "JB01",
        "202": "JB02",
        "203": "JB03",
        "204": "JB04",
        "205": "JB05",
        "206": "JB06",
        "207": "JB07",
        "208": "JB08",
        "209": "JB09",
        "210": "JB10",
        "211": "JB11",
        "212": "JB12",
        "213": "JB13",
        }
filename_convention = "{net}.{sta}.{loc}.{cha}.{typ}.{year}.{jday}"
outdir = "./"

# Assuming that all files in the directory are only waveforms
for fid in glob("*"):
    # NN.SSSSS.LL.CCC.T.YYYY.JJJ (SeisComP format)
    # network station location channel type year day-of-year
    net, sta, loc, cha, typ, year, jday = fid.split(".")

    try:
        new_net = network_dict[net]
        new_sta = station_dict[sta]
    except KeyError:
        print(f"{net} or {sta} not in dict")
        sys.exit(-1)

    # Change output name
    fid_out = filename_convention.format(net=new_net, sta=new_sta, loc=loc,
                                         cha=cha, typ=typ, year=year,
                                         jday=jday)
    print(f"{fid} -> {fid_out}")

    # Change actual stats
    st = read(fid)
    for tr in st:  # 1 tr per st expected
        tr.stats.network = new_net
        tr.stats.station = new_sta
        tr.stats.location = loc
        tr.stats.channel = cha
    
    st.write(f"{outdir}/{fid_out}", format="MSEED")


