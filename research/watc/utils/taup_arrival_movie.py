"""Watch a single phase traverse the globe"""
from obspy.taup import TauPyModel 
import matplotlib.pyplot as plt

phase_list = ["ttbasic"]
model = TauPyModel(model="prem")
for distance in range(0, 180, 1):
    arrivals = model.get_ray_paths(source_depth_in_km=0, 
                                   distance_in_degree=distance,
                                   phase_list=phase_list)
    arrivals.plot_rays(phase_list=phase_list, legend=True, show=False)
    plt.savefig(f"arrivals_{distance:0>3}.png")
    plt.close()
