"""
Plot RecordSections from ASDFDataSets
"""
import os
from glob import glob
from pyasdf import ASDFDataSet 
from pysep import RecordSection
from pysep.utils.io import read_asdfdataset


# for fid in ["/Users/chow/Work/work/akatom/nakversion/datasets/AK_A21K.h5"]:
for fid in ["/Users/chow/Work/work/akatom/nakversion/datasets/AK_D20K.h5"]:
# for fid in glob("/Users/chow/Work/work/akatom/nakversion/datasets/*.h5"):
    st, st_syn, windows = read_asdfdataset(fid, "i01s00")
    for comp in ["Z", "T"]:
        # Count windows
        n = 0
        for key in windows:
            if key.endswith(comp):
                n += len(windows[key])

        recsec = RecordSection(
                 st=st, 
                 #st_syn=st_syn, 
                 #windows=windows,
                 # min_period_s=30, 
                 # max_period_s=50, 
                 preprocess=None, 
                 scale_by="normalize", 
                 sort_by="distance", 
                 overwrite=True, 
                 components=comp, 
                 window_color="orange",
                 window_alpha=0.1,
                 linewidth=.7, 
                 tick_linewidth=0.,
                 title="",
                 y_axis_spacing=1,
                 xlim_s=[0, 400],
                 # move_out=4,
                 spine_top=False, 
                 spine_left=False, 
                 spine_right=False, 
                 show=False,
                 save=f"figures/{os.path.basename(fid).split('.')[0]}_{comp}_{n}.png"
                 ).run()    
