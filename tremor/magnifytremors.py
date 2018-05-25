"""indepth investigation of tremor detections by pyfreqscan
"""
import os
from obspy import UTCDateTime

def collect_files(date):
    date = UTCDateTime(date)
    fid_path = pathnames()['data'] + "TEROR/{y}/XX/*/pickle".format(y=date.year)
    files = glob.glob(os.path.join(fid_path,'*{}*'.format(date.julday)))
    if not files:
        print("No files found")
    return files
    
def 
    
                                                                  
    