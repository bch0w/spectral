"""
Creating waveform plots for illustrative purposes, so less labels,
bigger linesizes etc
"""
import os
import pyatoa
import pyasdf
from pyatoa import logger 
from pyflex import logger as pflogger

logger.setLevel("DEBUG")
pflogger.setLevel("DEBUG")

# General picks
picks = {"2014p952799": ["NZ.PXZ",   # forearc
                         "NZ.FWVZ",  # crustal
                         ],
         "2016p355601": ["NZ.BFZ"],  # crustal, 88km depth, no go
         "2016p356297": ["NZ.KHEZ",  # crustal through ocean to taranaki
                         "NZ.HIZ",   # crustal through land to taranaki
                         "NZ.PXZ"    # crustal wairarapa
                         ],
         "2017p012082": ["NZ.KHEZ",  # crustal through ocean poor fit
                         "NZ.TSZ"    # crustal through land poor fit
                         ],
         "2018p130600": ["NZ.TOZ",   # through TVZ
                         "NZ.MWZ",   # crustal north
                         ],
         "2019p738432": ["NZ.HIZ",   # through TVZ
                         "NZ.WIZ",   # crustal north
                         "NZ.MXZ",   # through forearc
                         "NZ.PUZ",   # through forearc
                         ],
         "3493233": ["NZ.BFZ"],      # from ocean
         "2016p275188": ["NZ.BKZ"],  # paper figure
         "3367989": ["NZ.WAZ",       # crossing horizontal
                     "NZ.VRZ",]
         }

# For paper1
picks = {"2014p952799": ["NZ.PXZ"],
         "2019p738432": ["NZ.WIZ"],
         "2017p012082": ["NZ.KHEZ"],
         "2018p130600": ["NZ.TOZ"]}

# picks = {"2019p738432": ["NZ.WIZ"]}  # MWZ also good
for min_period in [2, 10]:
    for event_id in picks:
        fid = os.path.join("./hiresdata", f"{event_id}.h5")
        assert(os.path.exists(fid)), fid

        config = pyatoa.Config(
            event_id=event_id,
            model_number="m00",
            min_period=min_period,
            max_period=30,
            filter_corners=4,
            component_list=["Z"],
            synthetics_only=False,
            rotate_to_rtz=False,
            unit_output="DISP",
            pyflex_preset="nznorth_10-30s",
            adj_src_type="cc",
            cfgpaths={"waveforms": [], "synthetics": './', "responses": []}
        )
        
        # Figure size based on which components were looking at
        if len(config.component_list) == 1:
            figsize = (10, 4)
        else:
            figsize = (10, 7)

        with pyasdf.ASDFDataSet(fid) as ds:
            for sta in picks[event_id]:
                station = f"{sta}.??.HH?"
                mgmt = pyatoa.Manager(config=config, ds=ds)
                mgmt.load(station_code=station, model=config.model_number)
                if len(config.component_list) == 1:
                    mgmt.st_obs = mgmt.st_obs.select(
                            component=config.component_list[0])
                mgmt.standardize()
                mgmt.preprocess()
                mgmt.window(fix_windows=False)
                mgmt.measure()

                mgmt.plot(show=False, 
                          save=(f"{config.event_id}_{sta}_"
                                f"{int(config.min_period)}-"
                                f"{int(config.max_period)}.png"), 
                          figsize=figsize, 
                          dpi=200, 
                          length_sec=200, 
                          fontsize=20,
                          axes_linewidth=3, 
                          linewidth=3., 
                          window_anno_fontsize=16, 
                          window_anno_height=0.0,
                          window_anno_rotation=0, 
                          window_color="dimgrey",
                          window_anno_fontcolor="k", 
                          window_anno_fontweight="roman", 
                          plot_waterlevel=False, 
                          plot_window_anno=True, 
                          plot_windows=True, 
                          plot_adjsrc=False, 
                          plot_stalta=False,
                          plot_arrivals=False,
                          legend=False)

