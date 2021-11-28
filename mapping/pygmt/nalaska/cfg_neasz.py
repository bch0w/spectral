# CONFIG FILE FOR THE NORTHEAST ALASKA SEISMIC ZONE (NEASZ)
# File IDs
output = "../figures/nalaska.png"
mt_fid = "/home/bchow/work/data/earthquakes/nalaska_mts_nolabel.txt"

# Toggle map features
moment_tensors = True

# Moment Tensor attributes
mt_scale = ".35c"  # Relative sizes of MT
convention = "mt"  # 'mt' for full moment tensor, 'dc' for closets double couple
mt_cmap = "seis"
mt_colorbar = True

# Colors
coast_color = "black"
bg_transparent = False

# Determine where we're mapping
region = [-151., -140., 66, 71.]
projection = "L-145/68/67/69/12"   # Lambert Conic Conformal Projection
