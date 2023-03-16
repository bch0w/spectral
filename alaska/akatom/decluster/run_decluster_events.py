"""
Using PySEP to curate the event catalog
"""
from obspy import read_events, read_inventory
from pysep import Declust, logger

logger.setLevel("DEBUG")
plot = False

cat = read_events("../metadata/nalaska_events.xml")
inv = read_inventory("../metadata/nalaska_stations.xml")

zedges = [0, 10, 35, 100]

declust = Declust(cat=cat, inv=inv)

# Kick out low-data events and threshold magnitudes by depth ranges
declust.threshold_catalog(zedges=zedges, min_mags=[3., 4., 5.], min_data=10)

# First pass: broad net to remove events with low data availability
dcl_cat = declust.decluster_events(choice="cartesian", zedges=zedges, 
                                   nkeep=[3, 3, 3], select_by="data_r", 
                                   nx=30, ny=30, plot=True)


# 2nd pass: more stringent net to get largest magnitudes only
dcl_cat2 = declust.decluster_events(cat=dcl_cat, choice="cartesian", 
                                    zedges=zedges, nkeep=[3, 2, 2], 
                                    select_by="magnitude", nx=15, ny=15,
                                    plot=False)


# Run srcrcv weight calculation
declust.calculate_srcrcv_weights(cat=dcl_cat, plot=True, show=True)

# Plot the original event catalog and the final event catalog
if plot:
    declust.plot(cat=cat, inv=inv, color_by="depth", vmin=0, vmax=30, 
                 cmap="inferno_r",
                 show=False, save="og_event_catalog.png", 
                 title=f"Original Event Catalog (N={len(cat)})")


    declust.plot(cat=cat, inv=inv, color_by="data", vmin=0, vmax=30, 
                 cmap="viridis",
                 show=False, save="og_data_avail.png", 
                 title=f"OG CAT DATA AVAIALABILITY (N={len(cat)})")
