STATIONS = 1
EVENTS = 0
LANDMARKS = 1
INTERFACE = 0
FAULTS = 0
BOUNDS = 0
TOPO = 1
RAYPATHS = 0
CITIES = 1

FIDS = {
    "EVENTS": "/Users/Chow/Documents/academic/vuw/data/events/"
              "decluster_60.xml",
    "STATIONS": "/Users/Chow/Documents/academic/vuw/data/specfem/stations/"
                "BEACON_GEONET",
    "INTERFACE": "/Users/Chow/Documents/academic/vuw/data/carto/interface/"
                 "williams_hikurangi_interface.npy",
    "FAULTS": "/Users/Chow/Documents/academic/vuw/data/carto/fault_coordinates/"
              "forest_mesh_gns_active_faults.txt",
    "TOPO": "/Users/Chow/Documents/academic/vuw/data/carto/topography/srtm30p/"
        "forest_mesh_topo_latlon.txt",
    "OUTPUT": "./figures/beacon_map.png"
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

FIGSIZE = (8, 10)
DPI = 150
AXIS_LINEWIDTH = 2.
DOMAIN_BOUNDS = {"lat_min": -42.5, "lat_max": -37.0,  
                 "lon_min": 173.0, "lon_max": 180}

DMN_KWARGS = {
    "facecolor": "None",
    "edgecolor": "k",
    "linewidth": 1.5,
    "linestyle": "--",
    "alpha": 1.,
    "zorder": 30,
}

MAP_KWARGS = {
    "map_corners": {"lat_min": -41.05, "lat_max": -38.95, 
                    "lon_min": 175.25, "lon_max": 177.6},
    "continent_color": "None",
    "lake_color": "None",
    "fill_color": "None",
    "coastline_linewidth": 3,
    "fontsize": 14,
    "area_thresh": 0,
    "scalebar": True,
    "scalebar_location": "lower-right",
    "scalebar_fontsize": 14.,
    "zorder": 20,
    "degrees": .5,
    "scalebar_length": 50
}

# STATION
STA_KWARGS = {
    "markersize": 75,
    "zorder": 100,
    "color": "w",
    "marker": "v",
    "edgecolor": "k",
    "annotate": True
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

