"""
Plotting Myanmar earthquake on HI stations
https://ds.iris.edu/wilber3/find_stations/11952284

starttime = UTCDateTime("2025-03-28T06:20:52")
"""
import os
from glob import glob
from obspy import read, UTCDateTime, Stream
from obspy.clients.fdsn import Client


ONE_HOUR = 60 * 60

comp = "Z"
sta_comp = "DESD"
stations = [f"HV.{sta_comp}.*.?H{comp}", f"IU.POHA.00.BH{comp}"]

path_hi = "/Users/chow/Data/nodes/2025-03_HAVO_EPSCoR"
path_out = "/Users/chow/Work/research/nasa_epscor24/HI_deployment_data/myanmar"

# Test
starttime = UTCDateTime("2025-04-26T14:20:54") - (ONE_HOUR * 0.5)
endtime = starttime + (ONE_HOUR * 0.5)

# Myanmar mww7.7 @ 10km depth
starttime = UTCDateTime("2025-03-28T06:20:52")
endtime = starttime + (ONE_HOUR + 4)

# Get permanent station data
c = Client("IRIS")
st = Stream()
for station in stations:
    net, sta, loc, cha = station.split(".")
    print(station)
    fid_st = os.path.join(path_out, f"{station}_data.ms")

    if not os.path.exists(fid_st):
        st_ = c.get_waveforms(network=net, station=sta, location=loc, 
                              channel=cha, starttime=starttime, 
                              endtime=endtime)
        st_.write(fid_st, "MSEED")
    else:
        st_ = read(fid_st)
    st += st_

# Get deployment data
for type_ in ["bare", "bucket", "windshield"]:
    subpath = f"{sta_comp}_{type_}"
    fullpath = os.path.join(path_hi, subpath)

    y = starttime.year
    m = f"{starttime.month:0>2}"
    d = f"{starttime.day:0>2}"
    for fid in glob(os.path.join(fullpath, f"*.*.{y}.{m}.{d}.*.{comp}.*")):
        print(fid)
        st_ = read(fid)
        st_[0].stats.network = "XX"
        st_[0].stats.station = sta_comp
        st_[0].stats.location = type_.upper()
        st += st_

st.trim(starttime, endtime)

st.detrend("linear")
st.detrend("simple")
st.filter("bandpass", freqmin=0.01, freqmax=10)

for tr in st:
    tr.write(os.path.join(path_out, f"{tr.get_id()}"), "MSEED") 

st.plot(equal_scale=False)





