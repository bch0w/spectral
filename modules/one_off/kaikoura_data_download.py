"""
Download data from GSN stations for Kaikoura Earthquake in SAC format
"""
from obspy import read_inventory
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import locations2degrees

# Get Kaikoura earthquake and information
c = Client("GEONET")
cat = c.get_events(eventid="2016p858000")
event = cat[0]
origintime = event.preferred_origin().time
endtime = origintime + 5000

# Get GSN station information
c = Client("IRIS")
inv = read_inventory("./GSNLIST.xml")
cha_code = ""
done_list = []
with open("source_receiver.txt", "w") as f:
	f.write("station,distance(deg),azimuth(deg),backazimuth(deg)\n")
	for net in inv:
		for sta in net:
			gcdist, az, baz = gps2dist_azimuth(
					lat1=event.preferred_origin().latitude,
					lon1=event.preferred_origin().longitude, lat2=sta.latitude,
					lon2=sta.longitude
					)
			dist_in_deg = locations2degrees(
					lat1=event.preferred_origin().latitude,
					long1=event.preferred_origin().longitude, lat2=sta.latitude,
					long2=sta.longitude
					)
			f.write("{sta},{dst:.2f},{az:.2f},{baz:.2f}\n".format(
				sta="{net}.{sta}".format(net=net.code, sta=sta.code),
				dst=dist_in_deg, az=az, baz=baz)
				)
			for cha in sta:
				if "{}.{}.{}".format(net.code, sta.code, cha.code) in done_list:
					continue
				try:
					st = c.get_waveforms(network=net.code, station=sta.code,
							location="*", channel=cha.code,
							starttime=origintime, endtime=endtime,
							attach_response=True
							)
					st.remove_response(output="DISP")
					st.trim(starttime=origintime, endtime=endtime)
					st.write(filename="{}.SAC".format(st[0].get_id()),
							format="SAC"
							)
					print("{}.{}.{} written".format(
						net.code, sta.code, cha.code)
						)
					done_list.append("{}.{}.{}".format(
						net.code, sta.code, cha.code)
						)
				except Exception as e:
					print("{}.{}.{}\n{}".format(
						net.code, sta.code, cha.code, e)
						)

