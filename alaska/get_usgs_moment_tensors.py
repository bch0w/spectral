"""
Get USGS moment tensors using ObsPy. Plot beachballs.
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from obspy import read_events, UTCDateTime, Catalog
from obspy.clients.fdsn import Client
from obspy.imaging.beachball import beach


def parse_event_id(event):
    """Just want the event tag and not the whole resource id"""
    return event.resource_id.id.split("eventid=")[1].split("&")[0]


# vvv PARAMETER SET
save_fid = "./usgs_alaska_cat.xml"  # Save the initial catalog
fm_save_fid = "./usgs_alaska_cat_w_fm.xml"  # Save the catalog with only focmecs
text_file = "./usgs_alaska_cat_w_fm.txt"  # Output text file with event info
plot = True  # Create a basic beachball plot using FM information
exclude_ak = False  # Don't consider any events tagged with AK (i.e., from AEC)

# Arguments to be sent to ObsPy Client.get_events()
kwargs = {
        "starttime": UTCDateTime("2000-01-01T00:00:00"),
        "endtime": UTCDateTime(),  # Present time
        "minlatitude": 56,
        "maxlatitude": 72,
        "minlongitude": 192,
        "maxlongitude": 222,
        "maxdepth": 200,
        "minmagnitude": 3,
        "includeallorigins": True,  # Required to get focal mechanism returns
        }

# Plotting parameters
cmap = plt.cm.jet_r
norm_a = .4  # controls relative size of earthquakes based on magnitude
norm_b = .75
xlim = [-168, -135]
ylim = [55, 72]
# ^^^ PARAMETER SET ^^^


# If this part has already been run, read existing catalog to save time 
if not os.path.exists(fm_save_fid):
    if not os.path.exists(save_fid):
        c = Client("USGS")
        cat = c.get_events(**kwargs)
        cat.write(filename=save_fid, format="QUAKEML")
    else:
        print(f"reading catalog from {save_fid}")
        cat = read_events(save_fid)

    print(f"Catalog has {len(cat)} events")

    # Not all events return focal mechanisms, kick those out
    events_with_fm = []
    for event in cat[:]:
        if event.focal_mechanisms and \
                   event.preferred_focal_mechanism().moment_tensor:
            if exclude_ak and ("ak" in parse_event_id(event)):
                continue
            events_with_fm.append(event)

    # Make the slimmed down catalog that only contains earthquakes with foc mecs
    cat = Catalog(events=events_with_fm)
    cat.write(filename=fm_save_fid, format="QUAKEML")
else:
    print(f"reading fm catalog from {fm_save_fid}")
    cat = read_events(fm_save_fid)
print(f"Catalog w/ focal mechanisms has {len(cat)} events")

# Plot moment tensors, colored by depth, scaled by magnitude
if plot:
    print("plotting crude beachball map")

    # Grab catalog information for relative colors and scaling
    xvals, yvals, depths, mags, event_ids = [], [], [], [], []
    for event in cat:
        xvals.append(event.preferred_origin().longitude)
        yvals.append(event.preferred_origin().latitude)
        depths.append(event.preferred_origin().depth * 1E-3)  # units km
        mags.append(event.preferred_magnitude().mag)
        event_ids.append(parse_event_id(event))
    print(f"{min(mags)} <= M <= {max(mags)}")
    print(f"{min(depths)} <= Z <= {max(depths)} km")

    # Might as well write out the above information to a text file
    arr = np.vstack((event_ids, xvals, yvals, depths, mags)).T
    with open(text_file, "w") as f:
        f.write("event_id lon lat depth(km) mag\n")
        np.savetxt(f, arr, fmt="%s")

    # Normalize magnitudes to get relative beachball sizes
    mags = np.array(mags)
    mags_norm = ((norm_b-norm_a) * (mags-mags.min()) / \
                                               (mags.max()-mags.min())) + norm_a

    # Create a colorscale based on depths
    max_depth = int(max(depths))
    normalize = mpl.colors.BoundaryNorm(range(0, max_depth, 5), cmap.N)

    f, ax = plt.subplots(1)
    failed_beachballs = 0
    for i, event in enumerate(cat[:]):
        try:
            fm = dict(event.preferred_focal_mechanism().moment_tensor.tensor)
            fm_list = [fm["m_rr"], fm["m_tt"], fm["m_pp"], 
                       fm["m_rt"], fm["m_rp"], fm["m_tp"]]
        except AttributeError:
            print(f"{parse_event_id(event)} failed plotting")
            failed_beachballs += 1

        x = event.preferred_origin().longitude
        y = event.preferred_origin().latitude
        z = event.preferred_origin().depth * 1E-3

        # Smaller events plotted on top
        b = beach(fm_list, xy=(x, y), width=mags_norm[i], 
                  facecolor=cmap(normalize(z)), linewidth=1, 
                  zorder=1/mags_norm[i])

        ax.add_collection(b)
  
    # Colorbar for depth range
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=normalize) 
    sm.set_array([])
    cbar = plt.colorbar(sm, shrink=0.45)
    cbar.set_label("depth [km]")

    # Text objects around the plot
    nevents = len(cat) - failed_beachballs
    plt.title(f"USGS Moment Tensor Catalog (N={nevents})\n"
              f"{kwargs['starttime'].year} <= YEAR <= {kwargs['endtime'].year}\n"
              f"{min(mags)} <= M <= {max(mags)}")
    plt.xlabel("Longitude [deg]")
    plt.ylabel("Latitude [deg]")

    plt.xlim(xlim)
    plt.ylim(ylim)

    ax.set_aspect(1)
    f.tight_layout()
    plt.show()


        




