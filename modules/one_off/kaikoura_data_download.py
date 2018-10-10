"""
Download data from GSN stations for Kaikoura Earthquake in SAC format
"""
from obspy import read_inventory, read
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import locations2degrees


def rotate_files_update_txt():
    """
    rotate all N/E to R/T and update the source_receiver.txt file in the process
    :return:
    """
    with open("source_receiver.txt","r") as f:
        lines = f.readlines()
    with open("source_receiver_new.txt","w") as f:
        f.write("station,distance(deg),azimuth(deg),backazimuth(deg)\n")
        for line in lines[1:]:
            stanet, dist, az, baz, flag = line.strip().split()
            north = glob.glob("{}.*.*N.SAC".format(stanet))
            east = glob.glob("{}.*.*E.SAC".format(stanet))
            if not north or not east:
                f.write("{stanet},{dst},{az},{baz},{flag}\n".format(
                    stanet=stanet, dst=dist, az=az, baz=baz,flag=1)
                )
            else:
                if len(north) > 1 or len(east) > 1:
                    print("lengths are longer than expected")
                    import ipdb; ipdb.set_trace()
                st = read(north[0])
                st += read(east[0])
                st.rotate(method="NE->RT", back_azimuth=float(baz))
                for tr in st:
                    tr.write(filename="{}.SAC".format(tr.get_id()),
                             format="SAC")
                f.write("{stanet},{dst},{az},{baz},{flag}\n".format(
                    stanet=stanet, dst=dist, az=az, baz=baz,flag=1)
                )



# Get Kaikoura earthquake and information
c = Client("GEONET")
cat = c.get_events(eventid="2016p858000")
event = cat[0]
origintime = event.preferred_origin().time
endtime = origintime + 5000

# Get GSN station information
with open("GSN_update.txt", "r") as f:
    lines = f.readlines()

c = Client("IRIS")
with open("source_receiver.txt", "a") as f:
    f.write("station,distance(deg),azimuth(deg),backazimuth(deg)\n")
    for line in lines[100:]:
        sta, net, lat, lon, _, _ = line.strip().split()
        gcdist, az, baz = gps2dist_azimuth(
                lat1=event.preferred_origin().latitude,
                lon1=event.preferred_origin().longitude, lat2=float(lat),
                lon2=float(lon)
                )
        dist_in_deg = locations2degrees(
                lat1=event.preferred_origin().latitude,
                long1=event.preferred_origin().longitude, lat2=float(lat),
                long2=float(lon)
                )
        f.write("{sta},{dst:.2f},{az:.2f},{baz:.2f},".format(
            sta="{net}.{sta}".format(net=net, sta=sta),
            dst=dist_in_deg, az=az, baz=baz)
            )
        try:
            st = c.get_waveforms(network=net, station=sta, location="*",
                    channel="BH?", starttime=origintime,
                    endtime=endtime, attach_response=True
                    )
            # check if more than 3 traces in stream
            if len(st) > 3:
                sampling_rates = []
                for tr in st:
                    print("{} {} {}".format(tr.stats.location, tr.stats.channel,
                        tr.stats.sampling_rate))
                    sampling_rates.append(tr.stats.sampling_rate)
                if 20 and 40 in set(sampling_rates):
                    st = st[3:]
                else:
                    import ipdb; ipdb.set_trace()

            # check if components in 1,2,Z, rotate
            complist = []
            for tr in st:
                if tr.stats.channel[-1] in ["1", "2"]:
                    inv = c.get_stations(network=st[0].stats.network,
                            station=st[0].stats.station,
                            location=st[0].stats.location,
                            channel="{}?".format(
                                tr.stats.channel[:2]),
                            starttime=origintime, endtime=endtime,
                            level="CHANNEL")
                    st.rotate(method="->ZNE", inventory=inv)
                    break

            st.remove_response(output="DISP")
            st.trim(starttime=origintime, endtime=endtime)
            for tr in st:
                tr.write(filename="{}.SAC".format(tr.get_id()), format="SAC")
            print("{} written".format(st[0].get_id()))
            f.write("1\n")
        except Exception as e:
            f.write("0\n")
            print("{}.{}\n{}".format(net, sta, e))

if __name__ == "__main__"