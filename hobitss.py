from obspy import read_inventory, read_events
from obspy.clients.fdsn import Client

# determine start and end times of project
inv = read_inventory('./datafiles/hobitss_stations.xml')
starts,ends = [],[]
for sta in inv[0]:
    starts.append(sta[0].start_date)
    ends.append(sta[0].end_date)
global_start = max(starts)
global_end = min(ends)

cat = read_events('./datafiles/hobitss_events.xml')

c = Client("IRIS")

for event in cat:
    eventid = str(event.resource_id).split('/')[1]
    starttime = event.origins[0].time - 100
    endtime = starttime + 600
    if (starttime <= global_start) or (endtime >= global_end):
        print(eventid,"skipped")
        continue
    print(event)
    EBS_data = c.get_waveforms(network="YH",
                               station="EBS*",
                               location="",
                               channel="*",
                               starttime=starttime,
                               endtime=endtime)
    LOBS_data = c.get_waveforms(network="YH",
                               station="LOBS*",
                               location="",
                               channel="HH*",
                               starttime=starttime,
                               endtime=endtime)
    EBS_data.plot(outfile='./output_plots/hobitss/{}_ebs.png'.format(eventid))
    LOBS_data.plot(outfile='./output_plots/hobitss/{}_lobs.png'.format(eventid))
