"""
Use TauPy to print out arrival times for a given source depth and station 
distance in degrees
"""
import sys
from obspy.taup import TauPyModel

phase_list = ["ttbasic"]
phase_list = ["Pdiff", "pPdiff", "PP"]

model = TauPyModel(model="prem")
arrivals = model.get_travel_times(source_depth_in_km=float(sys.argv[1]),
                                  distance_in_degree=float(sys.argv[2]),
                                  phase_list=phase_list
                                  )
print(arrivals)
arrivals = model.get_ray_paths(source_depth_in_km=float(sys.argv[1]),
                               distance_in_degree=float(sys.argv[2]),
                               phase_list=phase_list,)

arrivals.plot_rays()
