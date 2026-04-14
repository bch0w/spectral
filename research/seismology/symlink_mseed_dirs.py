import os
from glob import glob


dirs = ["/scale_wlg_nobackup/filesets/nobackup/gns03247/bchow/data/waveforms",
        "/home/chowbr/gns03247/project/bchow/data/waveforms"]

output = "/scale_wlg_nobackup/filesets/nobackup/gns03247/bchow/data/mseed"

for dir_ in dirs:
    for path in glob(os.path.join(dir_, "2*", "*")):
        year = os.path.basename(path) 

        for src in glob(os.path.join(path, "*")):
            path_out = os.path.join(output, year)
            if not os.path.exists(path_out):
                os.makedirs(path_out)

            dst = os.path.join(path_out, os.path.basename(src))          
 
            os.symlink(src, dst)
