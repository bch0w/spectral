"""
Use PyGMT to plot Northern Alaska domain with optional flags to fill up the 
space with random ish.

Notes:
	*) Pen options always follow the same format: width/color/texture
		- width goes from 0 (faint) to 12p (fattest)
		- color can be set with r/g/b, c/m/y/k, h-s-v, or by name
		- texture controls appearance, e.g., - or .. or ..-
"""
import os
import numpy as np
from pygmt import Figure, makecpt, config

overwrite_default_config = 0
# ============================= DEFAULT CONFIG =================================

# File IDs
output = "../figures/nalaska.png"
mt_fid = "/home/bchow/work/data/earthquakes/nalaska_mts_nolabel.txt"
landmarks_fid = "./landmarks.txt"

# Toggle map features
moment_tensors = True
inset = True
landmark_text = True

# Moment Tensor attributes
mt_scale = ".35c"  # Relative sizes of MT
convention = "mt"  # 'mt' for full moment tensor, 'dc' for closets double couple
mt_cmap = "seis"
mt_colorbar = True

# Colors
coast = "black"
land = "gray"
water = "white"
lakes = "gray"
bg_transparent = False

# Determine where we're mapping
region = [-168., -140., 64.5, 72.]  # Northern Alaska
frame = ["WSnE", "xa", "ya"]
projection = "L-155/68/67/69/12c"   # Lambert Conic Conformal Projection
map_scale = "g-144/71+w200"  # reference point and length

# Overwrite config parameters for specific look using .py files
if overwrite_default_config:
	print("OVERWRITING DEFAULT CONFIGURATION")
	from cfg_neasz import *

# ==============================================================================

# Set up standard map look
fig = Figure()
fig.coast(projection=projection, region=region, resolution="h",
          shorelines=[f"1/1.5p,{coast}", f"2/1p,{coast}"], 
          frame=frame, land=land, water=water, lakes=lakes,
		  area_thresh=10000)

# Plot moment tensors which should be in a file in the psmeca format from GCMT
# lon lat depth mrr mtt mpp mrt mrp mtp iexp name 
if moment_tensors and mt_fid:
	# Make the colormap to shade MTs by depth
	depths = np.loadtxt(mt_fid, usecols=2, dtype=float)
	makecpt(cmap=mt_cmap, series=[min(depths), max(depths), 1])
	
	mts = fig.meca(mt_fid, scale=mt_scale, convention="mt", component="dc",
				   C=True, verbose=True, L="1p/black")

	if mt_colorbar:
		fig.colorbar(position="n.95/.825+w-2c/.45c", frame='af+l"depth [km]"')

# Scalebar, need to use config to set the pen thickness of the scale
with config(MAP_TICK_PEN_PRIMARY=1.5):
	fig.basemap(region=region, projection=projection, #frame=["WSne", "gfa"],
				map_scale=map_scale)

# Landmark text
fig.text(textfiles=landmarks_fid)

# Plot the whole of Alaska as an inset to give context
if inset:
	with fig.inset(position="n0/.75+w5c/4c", margin=0):
		# Plot the whole of Alaska as the inset
		fig.coast(region=[-170, -135, 53, 72], projection="L-150/62/63/64/3.5c",
				  land="gray", water="white", shorelines="0p,black",
				  frame=False, resolution="h", area_thresh=100)

		# !!! THIS DOESNT WORK Put a box in the inset around Northern Alaska
		rectangle = [region[0], region[2], region[1], region[3]]
		fig.plot(data=rectangle, style="r+s", pen="1p,red")


# Plot the plate coupling / slip rate deficit with a no-label colorbar
# Save and show               
fig.savefig(output, transparent=bg_transparent)
fig.show(method="external")
