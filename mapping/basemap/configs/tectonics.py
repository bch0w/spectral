STATIONS = 0
EVENTS = 0
CITIES = 1
LANDMARKS = 1
INTERFACE = 0
FAULTS = 1
BOUNDS = 1
TOPO = 1
RAYPATHS = 0

FIDS = {
    "EVENTS": "",
    "STATIONS": "",
    "INTERFACE": "/Users/Chow/Documents/academic/vuw/data/carto/interface/"
                 "williams_hikurangi_interface.npy",
    "FAULTS": "/Users/Chow/Documents/academic/vuw/data/carto/fault_coordinates/"
              "forest_mesh_gns_active_faults.txt",
    "TOPO": "/Users/Chow/Documents/academic/vuw/data/carto/topography/srtm30p/"
            "forest_mesh_topo_latlon.txt",
    "OUTPUT": "./figures/tectonic_map.png"
}

# 1 for city, 0 for non-city
CITIES_DICT = {
    "Wellington": (-41.309, 174.759),
    r"Kaik$\mathrm{\bar{o}}$ura": (-42.405, 173.659),
    r"Taup$\mathrm{\bar{o}}$": (-38.6857, 176.0702),
    # "Hamilton": (-37.783, 175.252),
    "Auckland": (-36.903, 174.753),
    "Gisborne": (-38.665, 178.023, {"ha": "right"}),
    # "Te Araroa": (-37.636, 178.363, {"ha": "right"}),
    "East Cape": (-37.691, 178.539, {"ha": "right"}),
    "Whanganui": (-39.955, 175.028, {"ha": "right"}),
    r"P$\mathrm{\bar{o}}$rangahau": (-40.30, 176.61),
}

# Need to get newline working with raw string
mahia = r"M$\mathrm{\bar{a}}$hia" + "\nPeninsula"
LANDMARKS_DICT = {
    # "AUSTRALIAN\nPLATE": (-38.8, 172.5),
    # "PACIFIC PLATE": (-42., 177.5),
    # "ACCRETIONARY\nWEDGE": (-41.25, 176.3),
    # "CHATHAM RISE": (-42.85, 175),
    # "HIKURANGI\nPLATEAU": (-40.75, 178.2),
    "Mt. Taranaki": (-39.298, 174.063, {"ha": "center"}),
    "Bay of Plenty": (-37.5, 176.),
    "Cook Strait": (-41.6, 174.1),
    "Hawke Bay": (-39.6, 177.154),
    "Cook Strait": (-41.838, 174.525),
    "Mt. Ruapehu": (-39.283, 175.564, {"ha": "center"}),
    mahia: (-39.145, 177.9, {"ha": "right"}),
    "TVZ": (-38.33, 176)
}

FIGSIZE = (10, 12)
DPI = 150
AXIS_LINEWIDTH = 3.
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
    # "map_corners": {"lat_min": -42.75, "lat_max": -36.75,
    #                 "lon_min": 172.5, "lon_max": 178.8},
    "map_corners": {"lat_min": -42.9, "lat_max": -36.5,
                    "lon_min": 172, "lon_max": 179.25},
    "continent_color": "None",
    "lake_color": "None",
    "fill_color": "None",
    "coastline_linewidth": 2.,
    "fontsize": 14,
    "area_thresh": 10,
    "scalebar": True,
    "scalebar_location": "lower-right",
    "scalebar_fontsize": 14.,
    "zorder": 25,
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
    "textborder": "k",
    "linewidth": 2,
    "alpha": 0.5,
    "mark": False,
}

# CITIES
CTY_KWARGS = {
    "markersize": 50,
    "annotate": True,
    "fontsize": 11,
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
    "zorder": 99,
    "levels": [6, 9, 12, 15, 20, 30, 50, 100, 200, 300],
    "color": "k",
    "linestyle": "--",
    "alpha": .75,
    "linewidth": 2.5,
    "fontsize": 12,
    "fontweight": "heavy",
    "format_": "%.0f km",
}

# FAULTS
FLT_KWARGS = {
    "zorder": 98,
    "color": "k",
    "alpha": 1.,
    "linewidth": .7,
    "linestyle": "-"
}

TPO_KWARGS = {
    "zorder": 20,
    "cbar_shrink": .2,
    "land_segments": 200,
    "sea_segments": 200,
    "zero_col_val": .219
}
