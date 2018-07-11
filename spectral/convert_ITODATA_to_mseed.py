"""ito-san gave us sac data with no timestamps and in acceleration. set proper
time stamp, convert to velocity and save as mseed
"""
import glob
from obspy import read, UTCDateTime


timedict = {'2014p864702': UTCDateTime('2014-11-16T22:33:17'),
            '2014p051675': UTCDateTime('2014-01-20T02:52:45'),
            '2015p768477': UTCDateTime('2015-10-12T08:05:01'),
            '2014p715167': UTCDateTime('2014-09-22T14:41:22'),
            '2014p240655': UTCDateTime('2014-03-31T01:01:19'),
            '2016p859524': UTCDateTime('2016-11-14T00:34:22')
            }

sacfiles = glob.glob('*.s')
for sf in sacfiles:
    st = read(sf)
    event = sf.split('_')[2].split('.')[0]
    st[0].stats.starttime = timedict[event]
    newname = "{}_{}_vel.mseed".format(event,sf.split('_')[0])
    st.differentiate()
    for tr in st:
        tr.stats.station = sf.split('_')[0]
        tr.stats.network = "YH"
    st.write(newname,format="MSEED")
    
