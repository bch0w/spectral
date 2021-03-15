STATIONS = 1
EVENTS = 1
LANDMARKS = 0
INTERFACE = 1
FAULTS = 0
BOUNDS = 1

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

MAP_KWARGS = {
    "map_corners": {"lat_min": -42.75, "lat_max": -36.75, 
                    "lon_min": 172.5, "lon_max": 178.8},
    "continent_color": "whitesmoke",
    "lake_color": "azure",
    "fill_color": "azure",
    "coastline_linewidth": 2.,
    "fontsize": 12,
    "area_thresh": None,
    "scalebar": True,
    "scalebar_location": "lower-right",
    "scalebar_fontsize": 14.,
}

# STATION
STA_KWARGS = {
    "markersize": 90,
    "zorder": 100,
    "color": "w",
    "marker": "v",
    "edgecolor": "k"
}

# EVENT
EVT_KWARGS = {
    "cmap_name": "jet_r",
    "norm_a": 1.25E4,
    "norm_b": 2.75E4,
    "linewidth": 1.6,
    "zorder": STA_KWARGS["zorder"] - 2,
    "mag_scale": [4.5, 6.],
    "mag_legend": True,
    "cbar_shrink": 0.2,
    "cbar_fontsize": 15,
    "cbar_tickfontsize": 15,
    "cbar_labelpad": 17.5,
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
    "zorder": STA_KWARGS["zorder"] - 3,
    "levels": [0, 5, 10, 15, 25, 50, 100, 200],
    "color": "k",
    "linestyle": "-",
    "alpha": .75,
    "linewidth": 2.5,
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

