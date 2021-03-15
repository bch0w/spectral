STATIONS = 0
EVENTS = 0
LANDMARKS = 1
INTERFACE = 1
FAULTS = 1
BOUNDS = 1

FIDS = {
    "EVENTS": "",
    "STATIONS": "",
    "INTERFACE": "/Users/Chow/Documents/academic/vuw/data/carto/interface/"
                 "williams_hikurangi_interface.npy",
    "FAULTS": "/Users/Chow/Documents/academic/vuw/data/carto/fault_coordinates/"
              "forest_mesh_gns_active_faults.txt",
    "OUTPUT": "./figures/tectonic_map.png"
}

# 1 for city, 0 for non-city
LANDMARKS = {
    "Wellington": (-41.28664, 174.77557, 0),
    "Kaikoura": (-42.416665, 173.6833306, 0),
    "Taupo": (-38.6857, 176.0702, 0),
    "Mt. Taranaki": (-39.5, 174.0634, 1),
    "Bay of Plenty": (-37.5, 176., 1),
    "Cook Strait": (-41.6, 174.1, 1),
    "Mahia\nPeninsula": (-39.14, 177.90, 1),
    "Porangahau": (-40.30, 176.61, 0),
    "Australian Plate": (-38.4, 173.1, 1),
    "Pacific Plate": (-42., 177., 1)
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

