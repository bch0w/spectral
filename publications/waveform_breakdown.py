"""
Creating waveform plots for illustrative purposes, so less labels,
bigger linesizes etc
"""
import os
import pyatoa
import pyasdf

choice = ""

if choice == "inland":
    event_id = "2016p275188"
    sta_list = ["BKZ", "RTZ", "KHEZ", "KHZ"]
elif choice == "inland_2":
    event_id = "2019p754447"
    sta_list = ["HIZ"] # "TLZ", "TOZ"]
elif choice == "forearc_1":
    event_id = "2016p881118"
    sta_list = ["KNZ", "MXZ"]
elif choice == "forearc_2":
    event_id = "2015p768477"
    sta_list = ["KNZ", "PUZ", "MXZ"]
elif choice == "forearc_3":
    event_id = "3367989"
    sta_list = ["PUZ"] #"KNZ", "MXZ"]
elif choice == "forearc_back":
    event_id = "2019p738432"
    sta_list = ["PXZ"] # , "BFZ", "KHZ"]
elif choice == "volcano":
    event_id = "2014p240655"
    sta_list = ["TLZ"] # , "TOZ", "MKAZ"]
elif choice == "southern":
    event_id = "3620927"
    sta_list = ["NNZ", "MRZ", "KHZ"]
else:
    event_id = "2019p754447"
    sta_list = ["MRZ"]

config = pyatoa.Config(
    event_id=event_id,
    model_number="m00",
    min_period=10,
    max_period=30,
    filter_corners=4,
    # component_list=["Z"],
    synthetics_only=True,
    rotate_to_rtz=False,
    unit_output="DISP",
    pyflex_map="hikurangi",
    adj_src_type="cc_hikurangi",
    cfgpaths={"waveforms": [], "synthetics": './', "responses": []}
)

assert(os.path.exists(f"./hdf5/{config.event_id}.h5")), "event id doesn't exist"
ds = pyasdf.ASDFDataSet(f"./hdf5/{config.event_id}.h5")

for sta in sta_list:
    station = f"NZ.{sta}.??.HH?"
    mgmt = pyatoa.Manager(config=config, ds=ds)
    mgmt.populate(station_code=station, model=config.model_number)
    # mgmt.st_obs = mgmt.st_obs.select(component="Z")
    mgmt.preprocess()
    mgmt.run_pyflex()
    mgmt.run_pyadjoint()

    mgmt.plot_wav(show=False, save="clean2.png", figsize=(10,4), length_sec=180, 
                  fontsize=14`,
                  axes_linewidth=4, linewidth=2.5, window_anno_fontsize=0, 
                  window_anno_height=0.4, window_anno_rotation=0, 
                  window_anno_fontcolor="k", window_anno_fontweight="roman", 
                  plot_waterlevel=False, plot_window_anno=False, 
                  plot_windows=True, plot_adjsrc=True, plot_stalta=True,
                  legend=True)

