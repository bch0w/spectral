# File IDs
output = "figures/mahia_inset.png"
seis_fid = ("/Users/Chow/Documents/academic/vuw/tomo/interface/seismicity/"
            "mahia_eqs_geonet_2000_2021_z30_m2p5.csv_filtered.txt")

# Flags
plate_coupling = True
colorbar = False
interface = True
interface_bounds = False
sses = True
deep_sses = True
seamounts_ext = False
coast = True
trench = True
seismicity = True
seamounts = True
boxes = False
points = True

# Colors
coast_color = "white"
sse_color = "gold"
deep_sse_color = "greenyellow"
seamount_color = "limegreen"
seamount_ext_color = "lightslategray"

# Pens
sse_pen = f"3,{sse_color},solid"

# Fontsizes
sse_fontsize = "f12p"

# Determine where we're mapping
region = [177.42, 178.67, -39.75, -38.65]
projection = "U60S/12c"
interface_label = "l177/-39/178/-39.5"
map_scale = "g178.28/-39.65+c178.28/-39.65+w20"
