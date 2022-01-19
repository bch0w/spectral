STATIONS = 1
EVENTS = 1
LANDMARKS = 0
INTERFACE = 0
FAULTS = 0
BOUNDS = 0
TOPO = 0
RAYPATHS = 0
CITIES = 0

FIDS = {
    "EVENTS": "/Users/Chow/Documents/academic/vuw/south/prep/test_south.xml",
    # "EVENTS": "/Users/Chow/Documents/academic/vuw/south/prep/south_w_mt.xml",
    "STATIONS": "/Users/Chow/Documents/academic/vuw/south/prep/STATIONS_NZ_SOUTH",
    "INTERFACE": "/Users/Chow/Documents/academic/vuw/data/carto/interface/"
                 "williams_hikurangi_interface.npy",
    "FAULTS": "/Users/Chow/Documents/academic/vuw/data/carto/fault_coordinates/"
              "forest_mesh_gns_active_faults.txt",
    "TOPO": "/Users/Chow/Documents/academic/vuw/data/carto/topography/srtm30p/"
        "forest_mesh_topo_latlon.txt",
    "OUTPUT": "./figures/station_map.png"
}

# 1 for city, 0 for non-city
CITIES_DICT = {
        # "Porangahau": (-40.30, 176.61),
        # "Dannevirke": (-40.21, 176.09),
        # "Waipukurau": (-39.99, 176.55),
        # "Hastings": (-39.635, 176.831),
        # "Palmerston North": (-40.356, 175.613),
        # "Napier": (-39.51, 176.88),
}

LANDMARKS_DICT = {
        # "Castlepoint": (-40.90, 176.23),
        # "Cape Turnagain": (-40.49, 176.62),
        # "Cape Kidnappers": (-39.64, 177.095),


}

FIGSIZE = (12, 10)
DPI = 150
AXIS_LINEWIDTH = 2.
DOMAIN_BOUNDS = {"lat_min": -48., "lat_max": -40.,  
                 "lon_min": 165.0, "lon_max": 175}

DMN_KWARGS = {
    "facecolor": "None",
    "edgecolor": "k",
    "linewidth": 1.5,
    "linestyle": "--",
    "alpha": 1.,
    "zorder": 30,
}

MAP_KWARGS = {
    # Full map
    "map_corners": {"lat_min": -48.25, "lat_max": -39.9,  
                     "lon_min": 164.5, "lon_max": 175.8},
    "continent_color": "whitesmoke",
    "lake_color": "whitesmoke",
    "fill_color": "whitesmoke",
    "coastline_linewidth": 3,
    "fontsize": 14,
    "area_thresh": 0,
    "scalebar": True,
    "scalebar_location": "lower-right",
    "scalebar_fontsize": 14.,
    "zorder": 25,
    "degrees": 1,
    "scalebar_length": 100
}

# STATION
STA_KWARGS = {
    "markersize": 85,
    "zorder": 100,
    "color": "w",
    "marker": "v",
    "edgecolor": "k",
    "annotate": False
}

# EVENT
EVT_KWARGS = {
    "cmap_name": "jet_r",
    "norm_a": 2E4,
    "norm_b": 2.75E4,
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
    "markersize": 40,
    "annotate": True,
    "fontsize": 12,
    "zorder": 100,
    "color": "k",
    "edgecolor": "w",
    "fontcolor": "w",
    # "textborder": "k",
    "linewidth": 2,
    "alpha": 0.5,
    "mark": False,
}

# CITIES
CTY_KWARGS = {
    "markersize": 20,
    "annotate": False,
    "fontsize": 0,
    "fontcolor": "w",
    "textborder": "k",
    "fontweight": "normal",
    "zorder": 100,
    "color": "k",
    "linewidth": 2,
    "edgecolor": "w",
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

TPO_KWARGS = {
    "zorder": 20,
    "cbar_shrink": .2,
    "cmap_segments": 64,
    "zero_col_val": .219,
    "markersize": 1,
}

