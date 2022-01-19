STATIONS = 0
EVENTS = 0
RAYPATHS = 1
LANDMARKS = 0
INTERFACE = 0
FAULTS = 0
BOUNDS = 1
TOPO = 0
CITIES = 0

FIDS = {
    "EVENTS": "/Users/Chow/Documents/academic/vuw/data/events/"
              "decluster_60.xml",
    "STATIONS": "/Users/Chow/Documents/academic/vuw/data/specfem/stations/"
                "DECLUSTER",
    "INTERFACE": "/Users/Chow/Documents/academic/vuw/data/carto/interface/"
                 "williams_hikurangi_interface.npy",
    "FAULTS": "/Users/Chow/Documents/academic/vuw/data/carto/fault_coordinates/"
              "forest_mesh_gns_active_faults.txt",
    "OUTPUT": "./figures/raypath_map_b.png"
}

PAIRS = [
         ("2013p617227", "NZ.TOZ"),  # A
         ("2014p952799", "NZ.NTVZ"),  # B
         # ("2019p927023", "NZ.TMVZ"),
         # ("3367989", "NZ.NNZ"),
         # ("2013p614135", "NZ.MWZ"),
         # ("2016p105478", "NZ.ETVZ"), 
         ("2016p105478", "NZ.PUZ"),  # C
         # ("2016p105478", "NZ.WAZ"),
         # ("2016p356297", "NZ.KNZ"),
         ("2016p881118", "NZ.MWZ"),  # D
         # ("2019p927023", "NZ.NNZ"),
         ("2018p465580", "NZ.KHEZ"),  # E
         ("2019p738432", "NZ.KHZ"),  # F
         ("2019p754447", "NZ.HIZ"),  # G
         ("2019p927023", "NZ.VRZ"),]  # H

PAIRS = [PAIRS[4]]

# 1 for city, 0 for non-city
LANDMARKS = {
}

FIGSIZE = (4, 5)
DPI = 100
AXIS_LINEWIDTH = 2.
DOMAIN_BOUNDS = {"lat_min": -42.5, "lat_max": -37.0,  
                 "lon_min": 173.0, "lon_max": 178.5}

DMN_KWARGS = {
    "facecolor": "None",
    "edgecolor": "k",
    "linewidth": 1.5,
    "linestyle": "--",
    "alpha": 1.,
    "zorder": 30,
}

MAP_KWARGS = {
    "map_corners": {"lat_min": -42.75, "lat_max": -36.75, 
                    "lon_min": 172.5, "lon_max": 178.8},
    "continent_color": "w",
    "lake_color": "w",
    "fill_color": "w",
    "coastline_linewidth": 2.5,
    "fontsize": 14,
    "area_thresh": 0,
    "scalebar": True,
    "scalebar_location": "lower-right",
    "scalebar_fontsize": 14.,
    "zorder": 1
}

# STATION
STA_KWARGS = {
    "markersize": 60,
    "zorder": 100,
    "color": "w",
    "marker": "v",
    "edgecolor": "k"
}

# EVENT
EVT_KWARGS = {
    "cmap_name": "jet_r",
    "norm_a": 1.5E4,
    "norm_b": 2.5E4,
    "linewidth": 1.6,
    "zorder": STA_KWARGS["zorder"] - 2,
    "mag_scale": [4, 5., 6.],
    "mag_legend": True,
    "cbar_shrink": 0.2,
    "cbar_fontsize": 15,
    "cbar_tickfontsize": 15,
    "cbar_labelpad": 17.5,
    "cbar_linewidth": 2.,
}

# LANDMARKS
LMK_KWARGS = {
    "markersize": 100,
    "fontsize": 10,
    "zorder": STA_KWARGS["zorder"] - 1,
    "names": True,
    "color": "yellow",
    "alpha": 0.5,
    "marker": "o",
}

# PLATE INTERFACE
INT_KWARGS = {
    "zorder": STA_KWARGS["zorder"] - 6,
    "levels": [0, 5, 10, 15, 25, 50, 100, 200],
    "color": "k",
    "linestyle": "-",
    "alpha": 1.,
    "linewidth": 1.5,
    "fontsize": 12,
    "fontweight": "heavy",
    "format_": "%.0f km",
}

# FAULTS
FLT_KWARGS = {
    "zorder": STA_KWARGS["zorder"] - 5,
    "color": "k",
    "alpha": 1.,
    "linewidth": .5,
    "linestyle": ":"
}

