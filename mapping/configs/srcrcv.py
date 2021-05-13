STATIONS = 1
EVENTS = 1
LANDMARKS = 0
INTERFACE = 1
FAULTS = 0
BOUNDS = 1
TOPO = 0
RAYPATHS = 0

FIDS = {
    "EVENTS": "/Users/Chow/Documents/academic/vuw/data/events/"
              "decluster_60.xml",
    "STATIONS": "/Users/Chow/Documents/academic/vuw/data/specfem/stations/"
                "DECLUSTER",
    "INTERFACE": "/Users/Chow/Documents/academic/vuw/data/carto/interface/"
                 "williams_hikurangi_interface.npy",
    "FAULTS": "/Users/Chow/Documents/academic/vuw/data/carto/fault_coordinates/"
              "forest_mesh_gns_active_faults.txt",
    "OUTPUT": "./figures/srcrcv_map.png"
}

# 1 for city, 0 for non-city
LANDMARKS = {
}

FIGSIZE = (8, 10)
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
    "continent_color": "whitesmoke",
    "lake_color": "azure",
    "fill_color": "azure",
    "coastline_linewidth": 2.5,
    "fontsize": 14,
    "area_thresh": 100,
    "scalebar": True,
    "scalebar_location": "lower-right",
    "scalebar_fontsize": 14.,
    "zorder": 20
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
    # "levels": [3, 6, 9, 12, 15, 20, 30, 40, 50, 75, 100, 150, 200, 250 ,300],
    "levels": [6, 9, 12, 15, 20, 30, 50, 100, 200, 300],
    "color": "k",
    "linestyle": "--",
    "alpha": 1.,
    "linewidth": 1.,
    "fontsize": 1,
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

