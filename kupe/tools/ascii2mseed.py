"""Convert output ascii files from specfem into mseed files. Should be run
inside the folder you want to convert, i.e.

cd /path/to/mseeds/
python /path/to/ascii_to_mseed.py
"""
import os
import glob
import numpy as np
from obspy import UTCDateTime, read, Trace, Stream, read_events

dirname = os.getcwd()
all_files = glob.glob(os.path.join(dirname,'*.sem?'))
print("{} files found".format(len(all_files)))

cmtfilepath = glob.glob(os.path.join(dirname,'*CMTSOLUTION'))[0]
cmtfile = read_events(cmtfilepath,format='CMTSOLUTION')
if not cmtfile:
    sys.exit('CMTSOLUTION file required in folder')
starttime = cmtfile[0].origins[0].time


errors = 0
for fname in all_files:
    try:
        basename = os.path.basename(fname)
        output_name = "{}.mseed".format(basename)
        output_file = os.path.join(dirname,output_name)
        if os.path.exists(output_file):
            continue

        time = np.loadtxt(fname=fname, usecols=0)
        data = np.loadtxt(fname=fname, usecols=1)
        # assuming dt is constant after 3 decimal points
        delta = round(time[1]-time[0],3)

        network,station,channel,component = basename.split('.')
        stats = {"network":network,
                 "station":station,
                 "location":"",
                 "channel":channel,
                 "starttime":starttime,
                 "npts":len(data),
                 "delta":delta,
                 "mseed":{"dataquality":'D'}
                 }
        st = Stream([Trace(data=data,header=stats)])


        st.write(output_name,format="MSEED")
    except KeyboardInterrupt:
        sys.exit()
    except Exception as e:
        errors +=1
        print("Failure at ",e)

print("{} files converted".format(len(all_files)-errors))