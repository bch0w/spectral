"""
cgps data comes in full, trim down to the Chiapas earthquake
"""
import glob
from obspy import UTCDateTime


def read_txt(fid):
    with open(fid, 'r') as f:
        lines = f.readlines()
    return lines


def save_trimmed_text(fid, list_):
    new_fid = fid.split('.')[0] + "_trimmed." + fid.split('.')[1]
    with open(new_fid, 'w') as f:
        for entry in list_:
            f.write(entry)


cgps_files = glob.glob(
    "/Users/chowbr/Documents/subduction/spectral/"
    "tremor/gsnz18/raw/*_*.txt")
chiapas = UTCDateTime("2017-09-08T04:49:00.000Z")
plus_minus = 50
chiapas_minus = chiapas - plus_minus*(3600*24)
chiapas_plus = chiapas + (plus_minus+1)*(3600*24)

for fid in cgps_files:
    lines = read_txt(fid)
    day_clipped = [lines[0]]
    for i, line in enumerate(lines[1:]):
        if line[:10] == str(chiapas_minus)[:10]:
            index_start = i
        if line[:10] == str(chiapas_plus)[:10]:
            data = lines[index_start:i]
            save_trimmed_text(fid, data)
            break




