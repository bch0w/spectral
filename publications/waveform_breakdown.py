import pyatoa
import pyasdf

station = "NZ.MRZ.??.HH?"

config = pyatoa.Config(
    event_id="2019p754447",
    model_number="m00",
    min_period=10,
    max_period=30,
    filter_corners=4,
    rotate_to_rtz=False,
    unit_output="DISP",
    pyflex_map="hikurangi",
    adj_src_type="cc_hikurangi",
    cfgpaths={"waveforms": [], "synthetics": './', "responses": []}
)

ds = pyasdf.ASDFDataSet(f"{config.event_id}.h5")

mgmt = pyatoa.Manager(config=config, ds=ds)
mgmt.populate(station_code=station, model=config.model_number)
mgmt.preprocess()
mgmt.run_pyflex()
mgmt.run_pyadjoint()


mgmt.plot_wav(show=True, length_sec=200, axes_linewidth=4, linewidth=2.5, 
              fontsize=0, window_anno_fontsize=0, window_anno_height=0.4, 
              window_anno_rotation=0, window_anno_fontcolor="k",
              window_anno_fontweight="roman", legend=False)

