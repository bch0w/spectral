# Pulled from Prettyplot, need to clean up
path = ("/Users/chow/Work/research/gulkanaseis24/data/USGS_data/"
        "phelan_creek_stream_guage_2024-09-07_to_2024-09-14.txt")
assert(os.path.exists(path))

data = np.loadtxt(path, skiprows=28, usecols=[2,4], delimiter="\t", 
                    dtype=str)

# Convert units
if args.sg_units == "ft":
    height = height_ft
elif args.sg_units == "in":
    height = height_ft * 12
elif args.sg_units == "m":
    height = np.array([_ * 0.3048 for _ in height_ft.astype(float)])
elif args.sg_units == "cm":
    height = np.array([_ * 30.48 for _ in height_ft.astype(float)])
else:
    print("stream gage units `sg_units` should be in: ft, in, m, cm")
    sys.exit()

# Subset data where we are plotting waveforms to get the correct ylims
idx = np.where((times > xvals.min()) & (times < xvals.max()))

if args.sg_rel:
    height -= height[idx].min()

# Plot on the same axis as the waveform
twax = ax.twinx()

# # Plot the gradient, not very helpful
# dx = (times[1] - times[0]) * 86400
# grad_height = np.gradient(height[idx], dx)
# twax.plot(times[idx], grad_height, "o-", lw=1, c="C0", 
#         label="Phelan Cr. Gage", zorder=5, markersize=1.25, 
#         alpha=0.5)

twax.plot(times[idx], height[idx], "o-", lw=1, c="C0", 
            label="Phelan Cr. Gage", zorder=5, markersize=1.25, 
            alpha=0.5)

_ylabel = f"Stream Height [{args.sg_units}]"
if args.sg_rel:
    _ylabel = f"Relative {_ylabel}"
twax.set_ylabel(_ylabel, rotation=-90, labelpad=20)
