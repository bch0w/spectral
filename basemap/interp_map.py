import string
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyatoa.utils.srcrcv import lonlat_utm


def myround(x, base=5, choice='near'):
    """
    Round value x to nearest base, round 'up','down' or to 'near'est base

    :type x: float
    :param x: value to be rounded
    :type base: int
    :param base: nearest integer to be rounded to
    :type choice: str
    :param choice: method of rounding, 'up', 'down' or 'near'
    :rtype roundout: int
    :return: rounded value
    """
    if choice == 'near':
        roundout = int(base * round(float(x)/base))
    elif choice == 'down':
        roundout = int(base * np.floor(float(x)/base))
    elif choice == 'up':
        roundout = int(base * np.ceil(float(x)/base))

    return roundout


# Create a custom rotated line marker for a ruler effect on cross section
tn = mpl.markers.MarkerStyle(marker="|")
tn._transform = tn.get_transform().rotate_deg(-40)

tp = mpl.markers.MarkerStyle(marker="|")
tp._transform = tp.get_transform().rotate_deg(45)


# SIMULATION LIMITS
sim_x_min = 171312.
sim_x_max = 633468.
sim_y_min = 5286950.
sim_y_max = 5904080.

buff = 0
x_min = sim_x_min - buff
x_max = sim_x_max + buff
y_min = sim_y_min - buff
y_max = sim_y_max + buff

x_vals = np.arange(sim_x_min, sim_x_max, 1000)
y_vals = np.arange(sim_y_min, sim_y_max, 1000)

# COAST
coast = np.loadtxt("/Users/Chow/Documents/academic/vuw/data/carto/coastline/"
                   "coast_nznorth_utm60.txt")

# SLICES
# TRENCH_NORMAL = [-0.64, -0.76, 0.]
# TRENCH_PARALLEL = [-0.7, 0.7, 0.]
ORIGIN_NORMAL = [455763., 5547040., 0.]
LANDMARKS = {"Kaikoura": [226749., 5300535., 0.],
             "Wellington": [314007., 5426403., 0.,],
             "Flatpoint": [413092., 5433559., 0.,],
             "Castlepoint": [434724., 5471821., 0.,],
             "Akitio": [436980., 5511241., 0.],
             "Porangahau": [467051., 5538717., 0.,],
             "Elsthorpe": [484394., 5581561., 0.,],
             "Napier": [489374., 5626518., 0.,],
             "Mohaka": [507922., 5670909. ,0.,],
             "Mahia": [575567., 5665558., 0.,],
             "Gisborne": [588984., 5720001., 0.,],
             "Parallel": None,
             }

# Additional locations for reference
LOCATIONS = [(-38.6857, 176.0702, "Taupo"),
             (-39.2968, 174.0634, "Mt. Taranaki"),
             (-39.9333, 175.05, "Whanganui"),
             (-40.6463, 175.7054, "Ekatahuna"),
             (-38.1446, 176.2378, "Rotorua"),
             (-37.5139, 177.1811, "White Island"),
             (-39.2817, 175.5685, "Mt.Ruapehu"),
             (-37.7833, 175.2833, "Hamilton"),
             ]
for i, (name, origin) in enumerate(LANDMARKS.items()):
    plt.figure(figsize=(8,6), dpi=200)
    plt.scatter(coast[:,0], coast[:,1], s=.5, c='k')

    if (name == "Parallel"):
        b = (ORIGIN_NORMAL[1] - ORIGIN_NORMAL[0]) * .7
        y_vals = x_vals + b / .7
        plt.text(ORIGIN_NORMAL[0], ORIGIN_NORMAL[1], s=name, fontsize=8,
                 rotation=45, rotation_mode="anchor")
        t = tp
    else:
        # Create a marker for the origin point with annotation
        plt.scatter(origin[0], origin[1], c='y', s=40, marker='d', zorder=6)
        ab = string.ascii_uppercase[i]
        # plt.text(origin[0], origin[1], s=f"{ab}. {name}", fontsize=8,
        #          rotation=-40, rotation_mode="anchor")

        # Plot the cross-section line
        b = (.64 * origin[0] + .76 * origin[1])
        y_vals = (-.64 / .76) * x_vals + b / .76
        t = tn

    # Stop the line when it crosses the y-min val or the x-max val
    idy = np.where(y_vals < sim_y_min)[0]
    idx = np.where(y_vals > sim_y_max)[0]
    vals = np.unique(np.concatenate((idx, idy)))
   
    x_vals_ = np.delete(x_vals, vals)
    y_vals_ = np.delete(y_vals, vals)

    plt.plot(x_vals_, y_vals_, c="r", zorder=5)

    # Plot a scalebar along the lines
    dx_m = 50E3
    hypotenuse = np.sqrt((x_vals_ - x_vals_.min()) ** 2 +
                         (y_vals_ - y_vals_.max()) ** 2)
    max_length = np.sqrt((x_vals_.max() - x_vals_.min()) ** 2 +
                         (y_vals_.max() - y_vals_.min()) ** 2)

    for length in np.arange(0, myround(max_length, 1000, "up"), dx_m):
        val = min(hypotenuse, key=lambda x:abs(x-length))
        idx = np.where(hypotenuse == val)[0]
        plt.scatter(x_vals_[idx], y_vals_[idx], marker=t, c="r", s=80)

    for loc in LOCATIONS:
        lat, lon, name_ = loc
        x, y = lonlat_utm(lon, lat, -60)
        plt.scatter(x, y, c='orange', s=30, marker='v')
        plt.text(x, y, s=name_, fontsize=7)

    # Plot simulation domain
    plt.plot([sim_x_min, sim_x_max],[sim_y_min, sim_y_min], c='k', 
             linestyle='--')
    plt.plot([sim_x_min, sim_x_max],[sim_y_max, sim_y_max], c='k', 
             linestyle='--')
    plt.plot([sim_x_min, sim_x_min],[sim_y_min, sim_y_max], c='k', 
             linestyle='--')
    plt.plot([sim_x_max, sim_x_max],[sim_y_min, sim_y_max], c='k', 
             linestyle='--')

    # Bounds of the domain
    plt.xlim([x_min, x_max])
    plt.ylim([y_min, y_max])
    plt.gca().set_aspect(1)
    for axis in ["top", "bottom", "left", "right"]:
        plt.gca().spines[axis].set_linewidth(2)

    plt.gca().ticklabel_format(axis="both", style="sci", scilimits=(0,0))
    plt.title(f"{name} cross section (hashes every 50km)\n")
    plt.savefig(f"/Users/Chow/Documents/academic/vuw/forest/interpret/"
                f"refs/xsections/{ab}_{name.lower()}.png")
    plt.close()


