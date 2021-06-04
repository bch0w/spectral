"""
Use PyGMT to plot slip rate deficit

Things to add:
    1) offshore seamount locations
    2) plate convergence vector
    3) earthquake clusters
"""
import os
import numpy as np
from pygmt import Figure, makecpt
from xarray import DataArray

# Flags
plate_coupling = False
colorbar = False
interface = True
interface_bounds = False
sses = True
coast = True
trench = True
seismicity = True
points = True

# Determine where we're mapping
region = [173., 179, -42.5, -37.]
projection = "U60S/12c" 

# Set up standard map look
fig = Figure()
fig.coast(projection=projection, region=region, 
          shorelines=["1/1.5p,black", "2/1p,black"], frame=["WSne", "xa", "ya"])

# Skip because we don't want to see the grid lines
# fig.basemap(region=region, projection=projection, frame=["WSne", "gfa"])

# Set up colormap for plate coupling
cpt = "/Users/Chow/Documents/academic/vuw/data/carto/plate_coupling/rwb.cpt"
makecpt(cmap=cpt, series=[0, 1])

# Plot the plate coupling / slip rate deficit with a no-label colorbar
if plate_coupling:
    srd = "/Users/Chow/Documents/academic/vuw/data/carto/plate_coupling/hik.gmt"
    fig.plot(data=srd, cmap=True, close=True, color="+z")
    if colorbar:
        fig.colorbar(frame='af+l"Plate Coupling"')

# Plot the plate interface model as contours
if interface:
    intfc = "/Users/Chow/Documents/academic/vuw/data/carto/interface/williams_2013/supp_material_GMT/grid_exclude_wgs84.grd"
    contours = "/Users/Chow/Documents/academic/vuw/data/carto/interface/williams_2013/supp_material_GMT/contours.txt"
    fig.grdcontour(intfc, interval=contours, frame="f", # frame="a2f1g1", 
            annotation="40+f12", limit=[-100, 0], pen="thick,black,-", 
            label_placement="l173/-39/177/-42")

# Valid bounds for plate interface, not necessary
if interface_bounds:
    bounds = "/Users/Chow/Documents/academic/vuw/data/carto/interface/williams_2013/supp_material_GMT/bounds_wgs84.txt"
    fig.plot(data=bounds, columns=[0, 1], annotation="a2f1g1")

# Plot major SSEs as contours
if sses:
    sse = "/Users/Chow/Documents/academic/vuw/data/carto/sse_contours/filter_SSE2002_2014.grd"
    sse_color = "black"  # deeppink1
    fig.grdcontour(sse, pen=f"thicker,{sse_color},solid", limit=[100, 550], 
                   annotation=(100, "+f12"))

# Coastline towards the end so it plots over everything else
if coast:
    fig.coast(shorelines=["1/1.5p,black", "2/1p,black"], frame=["WSne", "xa", "ya"])

# Plot Hikurangi trench
if trench:
    trench_fid = "/Users/Chow/Documents/academic/vuw/data/carto/trench/hikurangi_trench_lonlat.txt"
    fig.plot(data=trench_fid, pen="thicker,black,solid")

if seismicity:
    base = "/Users/Chow/Documents/academic/vuw/publi/seamounts/figures/geosetting/"
    for fid in ["mahia_eqs_geonet_2000-2021_z20_m3.txt", "pora_eqs_geonet_2000-2021_z25_m3.txt"]:
        seis_fid = os.path.join(base, fid)
        fig.plot(data=seis_fid, style="x0.2c", color="black", pen="1.,black")

# Plot the locations of the seamounts as filled circles
if points:
    # Seamounts
    seamount_c = "limegreen"
    fig.plot(x=176.609, y=-40.300, style="c0.3c", color=seamount_c, pen="2.,black")  # Porangahau
    fig.plot(x=177.903, y=-39.134, style="c0.3c", color=seamount_c, pen="2.,black")  # Mahia Pen.

    # Other seamounts
    if False:
        seamount_c = "lightslategray"
        fig.plot(x=178.3, y=-40., style="c0.2c", color=seamount_c, pen="1.5,black")  # Rock Garden
        fig.plot(x=178.25, y=-41., style="c0.2c", color=seamount_c, pen="1.5,black")  # Bennett Knoll
        fig.plot(x=178.75, y=-40.03, style="c0.2c", color=seamount_c, pen="1.5,black")  # Offshore
        fig.plot(x=178.9, y=-39.1, style="c0.2c", color=seamount_c, pen="1.5,black")  # Bell SM1
        fig.plot(x=178.4, y=-39.6, style="c0.2c", color=seamount_c, pen="1.5,black")  # Bell 
        fig.plot(x=178.75, y=-38.8, style="c0.2c", color=seamount_c, pen="1.5,black")  # Bell 
        fig.plot(x=178.9, y=-38.4, style="c0.2c", color=seamount_c, pen="1.5,black")  # Bell 

    # Misc.
    # fig.plot(x=177.108, y=-40.515, style="d0.3c", color="cyan", pen="2.,black")  # Intraplate
    # fig.plot(x=177.792, y=-38.987, style="t0.35c", color="yellow", pen="2.,black")  # Morere
    # fig.plot(y=-40.393, x=176.905, style="+0.4c", pen="2.,black")  # Madden
    # fig.plot(y=-39.238, x=178.635, style="x0.4c", pen="2.,black")  # Poverty

# Save and show               
fig.savefig("figures/nznorth.png")
fig.show(method="external")
