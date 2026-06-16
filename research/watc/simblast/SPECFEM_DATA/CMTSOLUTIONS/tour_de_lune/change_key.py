import os
import sys
from glob import glob

key = sys.argv[1]
val = float(sys.argv[2])

for fid in glob("CMT*"):
    with open(fid) as f:
        lines = f.readlines()
    for i, line in enumerate(lines[:]):
        if line.startswith(key):
            lines[i] = f"{key}: {val}\n"
            break
    else:
        print(f"no key {key}")

    with open(fid, "w") as f:
        f.writelines(lines)

