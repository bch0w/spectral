"""30.1 written to overwrite all current stationxml files because the ones I had
didn't include all locations or components - wildcard in this function saves
the entire response file, obspy will figure out which one goes where later.
Probably a one time function
"""
import sys
import os
from obspy.clients.fdsn import Client
from obspy import UTCDateTime


station_list= ['BKZ','GKBS','GWTS','KFHS','MTHZ','MWZ','ORCS','ROTS','RTZ',
'TDHS','TPPS','TUDS','WAIS','WPWS','GISS','GTWS','HNPS','KNZ','MWDS','NSPS',
'PUZ','RPCS','TBAS','TIRS','TTHS','TUHS','WFSS']

net = 'NZ'
start = UTCDateTime('2015-01-01')
end = UTCDateTime('2016-01-01')
c = Client('GEONET')
for sta in station_list:
    print(sta)
    try:
        inv = c.get_stations(network=net,
                            station=sta,
                            location='*',
                            channel='*',
                            starttime=start,
                            endtime=end,
                            level='response')

        output_inv = '{net}_{sta}.xml'.format(net=net,sta=sta)
        filepath_inv = os.path.join(net,sta,output_inv)
        inv.write(filepath_inv,format="STATIONXML")
    except Exception as e:
        print(e)
        pass
