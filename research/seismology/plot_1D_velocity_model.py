"""
Read and plot the standard 1D velocity models from IRIS EMC. Allow plotting
multiple `choices` for comparison. Different models have different inherent 
colors to make it easier to visually separate.
Models are pulled directly from IRIS EMC links and will be downloaded
to your local machine
"""
import sys
import numpy as np
import matplotlib.pyplot as plt


choices = sys.argv[1:]
for choice in choices:
    assert choice in ["ak135f", "iasp91", "prem"]


# Choose which model to plot
depth_arrays = {}
plot_arrays = {}

# AK135f: depth, rho, vp, vs, qk, qm
if "ak135f" in choices:
    ak135f = "http://ds.iris.edu/files/products/emc/data/AK135F/AK135F_AVG.csv"
    depth_km, density, vp_kms, vs_kms, qk, qm = \
                        np.loadtxt(ak135f, delimiter=",").T

    depth_arrays["ak135f"] = depth_km
    plot_arrays["Vp_ak135f"] = vp_kms
    plot_arrays["Vs_ak135f"] = vs_kms
# IASP91: depth, radius, vp, vs
if "iasp91" in choices:
    iasp91 = "http://ds.iris.edu/files/products/emc/data/IASP91/IASP91.csv"
    depth_km, radius_km, vp_kms, vs_kms = np.loadtxt(iasp91, delimiter=",").T

    depth_arrays["iasp91"] = depth_km
    plot_arrays["Vp_iasp91"] = vp_kms
    plot_arrays["Vs_iasp91"] = vs_kms
# PREM Anisotropic: radius_km, depth_km, density, vpv, vph, vsv, vsh, eta, qm, qk
if "prem" in choices:
    prem_a = "http://ds.iris.edu/files/products/emc/data/PREM/PREM_1s.csv"
    radius_km, depth_km, density, vpv, vph, vsv, vsh, eta, qm, qk = \
            np.loadtxt(prem_a, delimiter=",").T

    depth_arrays["prem"] = depth_km
    plot_arrays["Vpv_prem"] = vpv
    plot_arrays["Vph_prem"] = vph
    plot_arrays["Vsv_prem"] = vsv
    plot_arrays["Vsh_prem"] = vsh

# Begin plotting
f, ax = plt.subplots(1, figsize=(5,8), dpi=125)
for name, arr in plot_arrays.items():
    if name.endswith("ak135f"):
        depth_km = depth_arrays["ak135f"]
        color_vp = "deepskyblue"
        color_vs = "lightcoral"
    elif name.endswith("iasp91"):
        depth_km = depth_arrays["iasp91"]
        color_vp = "springgreen"
        color_vs = "red"
    elif name.endswith("prem"):
        depth_km = depth_arrays["prem"]
        color_vp = "violet"
        color_vs = "coral"

    if name.startswith("Vp"):
        color = color_vp
    elif name.startswith("Vs"):
        color = color_vs
    
    plt.plot(arr, depth_km, color, label=name)

# Finalize plotting
plt.title(choices)
plt.xlabel("Velocity [km/s]")
plt.ylabel("Depth [km]")
plt.gca().invert_yaxis()
plt.grid()
plt.legend()
plt.ylim([depth_km.max(), -5])

plt.show()
