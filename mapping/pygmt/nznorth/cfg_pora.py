# File IDs
output = "figures/pora_inset.png"
seis_fid = ("/Users/Chow/Documents/academic/vuw/tomo/interface/seismicity/"
            "pora_eqs_geonet_2000_2021_z30_m2p5.csv_filtered.txt")
mt_fid = ("/Users/Chow/Documents/academic/vuw/tomo/interface/moment_tensors/"
          "mt_pora.csv")

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
seismicity = False
seamounts = True
boxes = False
points = True
moment_tensors = True

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
region = [176, 177.24, -40.9, -39.8]
projection = "U60S/12c"
interface_label = "l174.5/-39.25/178.5/-42" 
map_scale = "g176.9./-40.825+c176.9/-40.825+w20"
