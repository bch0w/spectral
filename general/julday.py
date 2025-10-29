"""CLI to get julday quickly"""
import sys
from obspy import UTCDateTime
time = UTCDateTime(sys.argv[1])
print(f"{time.year}-{time.julday:0>3}")
