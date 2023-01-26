"""
Using PySEP to curate the event catalog
"""
from obspy import read_events, read_inventory
from pysep import Declust, logger

logger.setLevel("DEBUG")


cat = read_events("../metadata/nalaska_events.xml")
inv = read_inventory("../metadata/nalaska_stations.xml")

zedges = [0, 10, 35, 100]

declust = Declust(cat=cat, inv=inv)

# Plot the original event catalog
declust.plot(cat=cat, inv=inv, color_by="depth", vmin=0, vmax=30, 
             cmap="inferno_r",
             show=False, save="og_event_catalog.png", 
             title=f"Original Event Catalog (N={len(cat)})")

declust.plot(cat=cat, inv=inv, color_by="data", vmin=0, vmax=30, 
             cmap="viridis",
             show=False, save="og_data_avail.png", 
             title=f"OG CAT DATA AVAIALABILITY (N={len(cat)})")

declust.threshold_events(zedges=zedges, min_mags=[3., 4., 5.], min_data=10)

# Select by data availability
newcat = declust.decluster_events(choice="cartesian", zedges=zedges, 
                                  nkeep=[3, 3, 3], select_by="data_r", 
                                  nx=30, ny=30, plot=True)

if False:
    polcat = declust.decluster_events(choice="polar", zedges=zedges, 
                                      nkeep=[3, 3, 3], select_by="data_r", 
                                      ntheta=32, plot=True)


declust.plot(cat=newcat, inv=inv, color_by="data", vmin=1, cmap="viridis",
             show=False, save="decluster_data_avail.png", 
             connect_data_avail=True,
             title="Declustered Event Catalog Data Availability")


