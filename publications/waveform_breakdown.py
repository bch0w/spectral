"""
Creating waveform plots for illustrative purposes, so less labels,
bigger linesizes etc
"""
import sys
import os
import pyatoa
import pyasdf
from pprint import pprint
from pyatoa import logger 
from pyflex import logger as pflogger
import matplotlib.pyplot as plt

logger.setLevel("INFO")
pflogger.setLevel("DEBUG")

min_period = int(sys.argv[1])
if min_period == 2:
    preset = "nznorth_2-30s"
else:
    preset = "nznorth_10-30s_plus"
show = True
comp = "Z"

# General picks
# For paper1 misfit gallery
# Geographic Variation
path = "./hiresdata"
picks = {
        "2017p012082": ["NZ.KHEZ", "A"],  # A
        "2018p130600": ["NZ.TOZ", "B"],   # B
        # "2019p738432": ["NZ.WIZ", "C"],   # C wrong
        "2014p952799": ["NZ.PXZ", "C"],   # D
        }
if min_period == 10:
    window_anno=("dT={cc_shift:.1f}s\n"
                 "dA={dlnA:.2f}\n"
                 "cc={max_cc:.2f}")
else:
    window_anno=("{cc_shift:3.1f}s\n"
                 "{dlnA:4.2f}\n"
                 "{max_cc:4.2f}")
window_anno_alternate=("{cc_shift:3.1f}s\n"
                       "{dlnA:4.2f}\n"
                       "{max_cc:4.2f}")


for event_id in picks:
    fid = os.path.join(path, f"{event_id}.h5")
    assert(os.path.exists(fid)), fid
    with pyasdf.ASDFDataSet(fid) as ds:
        config = pyatoa.Config(
            event_id=event_id,
            model="m00",
            min_period=min_period,
            max_period=30,
            filter_corners=4,
            component_list=[comp],
            synthetics_only=False,
            rotate_to_rtz=True,
            unit_output="DISP",
            pyflex_preset=preset,
            adj_src_type="cc",
            cfgpaths={"synthetics": './'}
        )
        
        # Figure size based on which components were looking at
        if len(config.component_list) == 1:
            figsize = (6.85, 2.5)
        else:
            figsize = (10, 7)

        sta = picks[event_id][0]
        label = picks[event_id][1]

        mgmt = pyatoa.Manager(config=config, ds=ds)
        mgmt.load(station_code=sta, path=config.model_number)
        mgmt.standardize()
        mgmt.preprocess()

        if len(config.component_list) == 1:
            mgmt.st_obs = mgmt.st_obs.select(
                    component=config.component_list[0])

        mgmt.window(fixed=False)

        # This is for the waveform gallery
        if event_id == "2014p952799" and sta == "NZ.PXZ":
            figsize = (6.25, 3.225)
            show_xaxis = True
        else:
            figsize = (6.25, 2.5)
            show_xaxis = False

        mp = mgmt.plot(show=False, 
                  save=False,
                  figsize=figsize, 
                  dpi=200, 
                  xlim_s=[15, 185], 
                  fontsize=20.,
                  axes_linewidth=2.5, 
                  linewidth=2.5, 
                  set_title=False,
                  percent_over=0.65,
                  window_anno=window_anno,
                  window_anno_alternate=window_anno_alternate,
                  window_anno_fontsize=16, 
                  window_anno_height=0.075,
                  alternate_anno_height=0.65,
                  window_anno_rotation=0, 
                  window_color="k",
                  window_anno_fontcolor="k", 
                  window_anno_fontweight="normal", 
                  plot_waterlevel=False, 
                  plot_window_anno=True, 
                  plot_windows=True, 
                  plot_adjsrc=False, 
                  plot_stalta=False,
                  plot_arrivals=False,
                  plot_rejected_windows=False,
                  plot_xaxis=show_xaxis,
                  plot_yaxis=False,
                  legend=False)
        # for win in mgmt.windows["Z"]:
        #     print(f"dT: {win.cc_shift * win.dt:.2f}\n"
        #             f"dA: {win.dlnA:.2f}\ncc: {win.max_cc_value:.2f}\n")
       
        lab_ = label
        if min_period == 2:
            lab_ += "'"

        mp.axes[0].text(x=0.02, y=0.825, s=f"{lab_}", 
                        verticalalignment="center", fontsize=30.,
                        transform=mp.axes[0].transAxes)
        mp.axes[0].tick_params(left=False, which='both')
        plt.savefig(f"./figures/{config.event_id}_{sta}_"
                    f"{int(config.min_period)}-"
                    f"{int(config.max_period)}_{comp}.png", figsize=figsize, 
                    dpi=200)
        # test = input("waiting")
