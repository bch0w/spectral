"""
Symlink FORCESOLUTION files into new sub-directories to avoid having repeat 
files. This also allows for single changes to the source files to propagate
to the other directories, e.g., in the case when one of the my sources had 
a weird float error due to its location
"""
import os
from glob import glob


for fid in glob("sources_*.txt"):
    lines = open(fid).readlines()
    name = os.path.splitext(fid)[0][8:].upper()
    dir_ = f"FORCESOLUTIONS_{name}"
    os.makedirs(dir_)
    for line in lines:
        if line.startswith("#"):
            continue
        name = line.strip()
        src = f"../FORCESOLUTIONS_ALL/{name}"
        dst = f"./{dir_}/{name}"
        os.symlink(src, dst)


