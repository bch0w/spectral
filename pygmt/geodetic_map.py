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


# ============================= DEFAULT CONFIG =================================
# File IDs
output = "figures/geodetic_revised.png"
# output = "/Users/Chow/Documents/academic/vuw/publi/seamounts/figures/geosetting/figures/geodetic_only.png"
seis_fid = "/Users/Chow/Documents/academic/vuw/data/events/geonet_forest_north_mgt6/earthquakes_norm.txt"
mt_fid = None

# Flags
plate_coupling = True
colorbar = False
interface = True
interface_bounds = False
sses = True
deep_sses = True
seamounts_ext = False
tsunami_eq = True
coast = True
trench = True
seismicity = True
seamounts = True
boxes = False
points = True
moment_tensors = True

# Colors
coast_color = "black"
sse_color = "gold"
tsunami_eq_color = "gold"
deep_sse_color = "slategray1"  # "deeppink1"
seamount_color = "limegreen"
seamount_ext_color = "lightslategray"
box_color = "yellow"
bg_transparent = False

# Pens
sse_pen = f"1.5,{sse_color},solid" 
deep_sse_pen = f"1.5,{deep_sse_color},solid"

# Fontsizes 
sse_fontsize = "f8p"

# Determine where we're mapping
region = [173., 179, -42.5, -37.]
projection = "U60S/12c" 
interface_label = "l173/-39/177/-42"
map_scale = "g178/-42+c178/-42+w100"

# Overwrite config parameters for specific look using .py files, comment to skip
# from mahia_cfg import *
from cartoon_cfg import *
# ==============================================================================


# Set up standard map look
fig = Figure()
fig.coast(projection=projection, region=region, 
          shorelines=[f"1/1.5p,{coast_color}", f"2/1p,{coast_color}"], 
          frame=["WSne", "xa", "ya"])
          # land="white", water="white")

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
    intfc = ("/Users/Chow/Documents/academic/vuw/data/carto/interface/"
             "williams_2013/supp_material_GMT/grid_exclude_wgs84.grd")
    contours = ("/Users/Chow/Documents/academic/vuw/data/carto/interface/"
                "williams_2013/supp_material_GMT/contours.txt")
    fig.grdcontour(intfc, interval=contours, frame="f", # frame="a2f1g1", 
                   annotation="40+f12", limit=[-100, 0], pen="thick,black,-", 
                   label_placement=interface_label
                   )

# Valid bounds for plate interface, not necessary
if interface_bounds:
    bounds = ("/Users/Chow/Documents/academic/vuw/data/carto/interface/"
              "williams_2013/supp_material_GMT/bounds_wgs84.txt")
    fig.plot(data=bounds, columns=[0, 1], annotation="a2f1g1")

# Coastline towards the end so it plots over everything else
if coast:
    fig.coast(shorelines=["1/1.5p,black", "2/1p,black"], 
              frame=["WSne", "xa", "ya"])

# Plot Hikurangi trench
if trench:
    trench_fid = ("/Users/Chow/Documents/academic/vuw/data/carto/trench/"
                  "hikurangi_trench_lonlat.txt")
    fig.plot(data=trench_fid, pen="thicker,black,solid",
             style="f1c/0.3c+r+t")  # carrots denoting subduction direction


# Plot major SSEs as contours
if sses:
    sse = ("/Users/Chow/Documents/academic/vuw/data/carto/slow_slip_events/"
           "filter_SSE2002_2014.grd")
    fig.grdcontour(sse, pen=sse_pen, limit=[100, 550],
                   annotation=(100, sse_fontsize))
# Deep 2006, 2008 slip
if deep_sses:
    deep = ("/Users/Chow/Documents/academic/vuw/data/carto/slow_slip_events/"
            "sum_06_08_deep.txt")
    fig.contour(data=deep, levels=20, pen=deep_sse_pen, annotation=(20, "+f8"))

if seismicity:
    fig.plot(data=seis_fid, style="a", color="white", pen="1.,black")  # style='c'

if moment_tensors and mt_fid:
    mts = open(mt_fid, "r").readlines()
    for mt in mts[1:]:
        _, _, lat, lon, strike, dip, rake, *_ = mt.split(",")
        *_, ml, mw, m0 = mt.split(",")
        fig.meca(dict(strike=float(strike), dip=float(dip), rake=float(rake), 
                      magnitude=float(mw)/2),
                 scale="1c", longitude=float(lon), latitude=float(lat),
                 depth=0
                 )

# Plot the locations of the 1947 Tsunami earthquakes
if tsunami_eq:
    # M 7.0 Gisborne, Mar 26 1947 (Offshore Poverty Bay)
    # Markersize taken from normalized earthquake list
    fig.plot(x=178.8, y=-38.85, style="a0.4099c", color="mediumpurple1", 
             pen="1,black") 
    # M 6.9 Gisborne, May 17 1947 (Offshore Tolaga Bay)
    fig.plot(x=178.87, y=-38.42, style="a0.3989c", color="deeppink", 
             pen="1,black") 

# Draw rectangles related to insets
if boxes:
    fig.plot(data=np.array([[176.,-40.9,177.25,-39.8]]), style='r+s', 
             pen=f"2p,{box_color},-")
    fig.plot(data=np.array([[177.42,-39.75,178.67,-38.65]]), style='r+s', 
             pen=f"2p,{box_color},-")

# No frame because we don't want to see the grid lines
# map_scale > g=ref point, c=set scale at, w=length
fig.basemap(region=region, projection=projection, #frame=["WSne", "gfa"],
            map_scale=map_scale)

# Plot the locations of the seamounts 
if seamounts:
    # Porangahau
    fig.plot(x=176.609, y=-40.300, style="x0.6c", color=seamount_color, 
             pen="4.5,black") 
    fig.plot(x=176.609, y=-40.300, style="x0.5c", color=seamount_color, 
             pen="2.5,lawngreen") 
    # Mahia Pen.
    fig.plot(x=177.903, y=-39.134, style="x0.6c", color=seamount_color, 
            pen="4.5,black") 
    fig.plot(x=177.903, y=-39.134, style="x0.5c", color=seamount_color, 
             pen="2.5,lawngreen") 

# Other seamounts
if seamounts_ext:
    style = "c0.2c"
    fig.plot(x=178.3, y=-40., style=style, color=seamount_ext_color, 
             pen="1.5,black")  # Rock Garden
    fig.plot(x=178.25, y=-41., style=style, color=seamount_ext_color,
             pen="1.5,black")  # Bennett Knoll
    fig.plot(x=178.75, y=-40.03, style=style, color=seamount_ext_color,
             pen="1.5,black")  # Offshore
    fig.plot(x=178.9, y=-39.1, style=style, color=seamount_ext_color, 
             pen="1.5,black")  # Bell SM1
    # fig.plot(x=178.4, y=-39.6, style=style, color=seamount_ext_color, 
    #          pen="1.5,black")  # Bell 
    fig.plot(x=178.75, y=-38.8, style=style, color=seamount_ext_color, 
             pen="1.5,black")  # Bell 
    fig.plot(x=178.9, y=-38.4, style=style, color=seamount_ext_color, 
             pen="1.5,black")  # Bell 

if points:
    fig.plot(x=177.108, y=-40.515, style="x0.5c", color="cyan", 
             pen="3.5,black")  # Intraplate
    fig.plot(x=177.108, y=-40.515, style="x0.4c", color="cyan", 
             pen="2.,cyan")  # Intraplate
    fig.plot(x=177.792, y=-38.987, style="+0.5c", color="yellow", 
             pen="3.5,black")  # Morere
    fig.plot(x=177.792, y=-38.987, style="+0.4c", color="yellow", 
             pen="2.,yellow")  # Morere
    fig.plot(y=-40.393, x=176.905, style="c0.01c", pen="1.,black")  # Madden
    fig.plot(y=-39.238, x=178.635, style="c0.01c", pen="1.,black")  # Poverty

# Save and show               
fig.savefig(output, transparent=bg_transparent)
fig.show(method="external")
