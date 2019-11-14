"""
A short Python script used to download GeoNet data

This requires ObsPy (obspy.org) and Python3

Information:
1) different FDSN clients:
https://docs.obspy.org/packages/obspy.clients.fdsn.html
2) filtering options:
https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.filter.html?highlight=filter#obspy.core.stream.Stream.filter
3) available file formats:  
https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.write.html#obspy.core.stream.Stream.write
"""
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# vvv  SET PARAMETERS HERE vvv
# Station
client = "GEONET"  # GeoNet if you want NZ data
network = "NZ"  # this is usually static
station = "BFZ"   
location = "*"  # GeoNet broadband locations are usually '10'
channel = "HH?"  # wildcard component to get 3 components

# Event
# starttime = UTCDateTime("2014-01-20T02:52:45Z")  # e.g. Eketahuna M6.2
# endtime = starttime + 90  # end of waveform, X seconds are starttime
starttime = UTCDateTime("2016-11-13T11:02:56Z")  # e.g. Kaikoura M7.8
endtime = starttime + 500  # end of waveform, X seconds are starttime

# Processing
remove_response = True  # if you want to remove the instrument response
min_period = 10  # for bandpass filter, if 'None', no filter applied
max_period = 30

# Output
output_format = "MSEED"  # also SAC, SEGY, SU etc.
# ^^^ SET PARAMETERS HERE ^^^


# Set the filename for saving the waveforms
fid_out = "{net}_{sta}_{year}_{jday}".format(
        net=network, sta=station, year=starttime.year, jday=starttime.julday)

# Get waveform data as an Obspy Stream object
c = Client(client)
st = c.get_waveforms(network=network, station=station, location=location,
                     channel=channel, starttime=starttime, endtime=endtime,
                     attach_response=True)

# Write the raw data
st.write(fid_out + "_raw.{}".format(output_format), format=output_format)

# Remove response
if remove_response:
    st.remove_response()

# Preprocessing functionality
st.detrend("linear")  # detrend the data to remove any very-long period signal
st.detrend("demean")  # from the data
st.taper(max_percentage = 0.05)  # taper the ends since we cut the data
# Filter the data 
if min_period:
    st.filter("bandpass", freqmin=1/max_period, freqmax=1/min_period)
st.detrend("linear")  # detrend and taper again incase filtering created any
st.detrend("demean")  # spurious signals
st.taper(max_percentage = 0.05)

# print the stream object for information
print(st)

# plot the stream and save it to the current directory
st.plot(outfile="./{}.png".format(fid_out))

# Save the data based on the User-defined file format
st.write(fid_out + "_processed.{}".format(output_format), format=output_format)





