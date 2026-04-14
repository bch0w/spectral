"""
Trying to make sure PySEP is returning the correct amplitdues. Compare
figure generated here against:

https://service.iris.edu/irisws/timeseries/1/query?net=XV&sta=F2TN&cha=HHZ&start=2018-10-03T03:29:00&end=2018-10-03T03:31:00&scale=AUTO&lp=1&format=plot&loc=--

"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

start = UTCDateTime("2018-10-03T03:26:17.550000Z")
end = UTCDateTime("2018-10-03T03:31:37.540000Z")

c = Client("IRIS")
st = c.get_waveforms("XV", "F2TN", "*", "HH?", start, end, attach_response=True)

pre_filt = [.00232, .00546, 25, 50]
st.remove_response(pre_filt=pre_filt, output="VEL", water_level=10000)
st.filter("lowpass", freq=1)

for tr in st:
    print(f"{tr.get_id()}: {tr.data.max():.2E}")

st.plot(outfile="2018-10-03_XV_F2TN.png", equal_scale=False,
        starttime=UTCDateTime("2018-10-03T03:29:00.000000Z"),
        endtime=UTCDateTime("2018-10-03T03:31:00.000000Z"))

