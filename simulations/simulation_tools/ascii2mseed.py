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
all_files = glob.glob(os.path.join(dirname,'*.sem?'))+\
									glob.glob(os.path.join(dirname,'*.adj'))
print("{} files found".format(len(all_files)))

cmtfilepath = glob.glob(os.path.join(dirname,'*CMTSOLUTION'))
if cmtfilepath:
    cmtfile = read_events(cmtfilepath[0],format='CMTSOLUTION')
    starttime = cmtfile[0].origins[0].time
else:
    print("ACHTUNG: No CMTSOLUTION file; starttime set as 2000-01-01T00:00:00")
    starttime = UTCDateTime('2000-01-01T00:00:00')
    cont = input("Continue? (y/[n])")
    if cont == 'n':
        import sys
        sys.exit()
    

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

        try:
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
        except ValueError:
            print('nonstandard filename - cannot write header info')
            st = Stream([Trace(data=data)])

        st.write(output_name,format="MSEED")
    except KeyboardInterrupt:
        sys.exit()
    except Exception as e:
        errors +=1
        print("Failure at ",e)

print("{} files converted".format(len(all_files)-errors))

deletecheck = input('Delete ascii files? (y/[n])')
if deletecheck == 'y':
	for fname in all_files:
		os.remove(fname)
