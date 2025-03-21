import os
import matplotlib.pyplot as plt
from obspy import UTCDateTime, read, read_inventory
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth, locations2degrees
from obspy.taup import TauPyModel


# 2024 Noto
origintimes = {"noto": UTCDateTime("2024-01-01T07:10:00"),
               "hualien": UTCDateTime("2024-04-03T23:58:11"),
               "kaikoura": UTCDateTime("2016-11-13T11:02:00")
               }

locations = {"noto": (37.488, 137.271), 
             "hualien": (23.819, 121.562),
             "kaikoura": (-42.737, 173.054)
             } 

stations = {"noto": ["IU.MAJO.*.BH?", "II.ERM.*.BH?"],
            "hualien": ["IU.TATO.*.BHZ", "HK.HKPS.BH?"],
            "kaikoura": ["NZ.WEL.20.HN?", "NZ.WEL.10.HH?"]  # 20=strong motion
            }

clients = {"noto": "IRIS", "hualien": "IRIS", "kaikoura": "GEONET"}

choice = "kaikoura"
duration = 90.
reset = True

c = Client(clients[choice])
origintime = origintimes[choice]
for station in stations[choice]:
    net, sta, loc, cha = station.split(".")

    fid = f"{choice}_{station}"
    inv_fid = f"inv_{fid}.xml"
    st_fid = f"st_{fid}.ms"

    if not reset and os.path.exists(inv_fid):
        inv = read_inventory(inv_fid)
    else:
        inv = c.get_stations(starttime=origintime, endtime=origintime+duration,
                             network=net, station=sta, location=loc, 
                             channel=cha, level="response")
        inv.write(inv_fid, format="STATIONXML")

    if not reset and os.path.exists(st_fid):
        st = read(st_fid)
    else:
        st = c.get_waveforms(network=net, station=sta, location=loc, 
                             channel=cha, starttime=origintime, 
                             endtime=origintime + duration
                             )
        st.write(st_fid, format="MSEED")
    
    eq_lat, eq_lon = locations[choice]
    rcv_lat = inv[0][0].latitude   # assuming inv is per-station
    rcv_lon = inv[0][0].longitude   # assuming inv is per-station

    dist_m, az, baz = gps2dist_azimuth(lat1=eq_lat, lon1=eq_lon,
                                       lat2=rcv_lat, lon2=rcv_lon)

    start = st[0].stats.starttime
    end = st[0].stats.endtime
    for tr in st:
        if tr.stats.starttime > start:
            start = tr.stats.starttime
        if tr.stats.endtime < end:
            end = tr.stats.endtime
    st.trim(start, end)


    st.remove_response(output="VEL", inventory=inv)
    st.rotate(method="->ZNE", inventory=inv)
    st.rotate(method="NE->RT", back_azimuth=baz)
    st.filter("bandpass", freqmin=1, freqmax=100)

    # TauP arrivals
    model = TauPyModel(model="iasp91")
    dist_deg = locations2degrees(lat1=eq_lat, long1=eq_lon, 
                                 lat2=rcv_lat, long2=rcv_lon)
    arrivals = model.get_travel_times(source_depth_in_km=15, 
                                      distance_in_degree=dist_deg,
                                      phase_list=["P", "S"])

    
    f, axs = plt.subplots(3, sharex=True, dpi=100)
    for i, component in enumerate(["Z", "R", "T"]):
        tr = st.select(component=component)[0]
        axs[i].plot(tr.times(), tr.data, c="k", lw=0.5, label=tr.get_id())
        plotted = []
        C = 0
        for arrival in arrivals:
            if arrival.name not in plotted:
                axs[i].axvline(arrival.time, c=f"C{C}")
                plotted.append(arrival.name)
                C += 1
        axs[i].legend()

    axs[2].set_xlabel("Time [s]")
    axs[1].set_ylabel("Amplitude")

    plt.tight_layout()
    plt.savefig(f"{fid}.png")
    plt.show()



