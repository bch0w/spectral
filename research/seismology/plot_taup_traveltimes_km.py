"""
Plot TauP travel-time curves with distance in kilometers rather than
degrees, focused on the near-source distance range where shallow
crustal discontinuities produce triplications.

Mirrors :func:`obspy.taup.plot_travel_times` but reads the precomputed
ray branches directly so triplication branches are not joined across
shadow zones, and converts distance to kilometers using the model's
own planet radius.
"""
import matplotlib.pyplot as plt
import numpy as np
from obspy.taup import TauPyModel
from obspy.taup.seismic_phase import SeismicPhase
from obspy.taup.utils import parse_phase_list

DEFAULT_MODEL = (
    "/Users/prof/Repos/spectral/research/seismology/taup_models/"
    "ak135f_upper_crust.npz"
)
# "P"/"S" are omitted: within the near-source crustal window this
# model is built for, they retrace the same branches as "Pg"/"Sg" and
# so would just draw an identical, fully overlapping line.
DEFAULT_PHASES = ["Pg", "Pn", "PmP", "Sg", "Sn", "SmS"]

# Each P-family phase gets its own cool sequential colormap, each
# S-family phase its own warm one; a phase's triplication legs are
# then different shades sampled from that one colormap.
COOL_CMAPS = ["Blues", "Greens", "Purples", "BuGn", "GnBu", "PuBu"]
WARM_CMAPS = ["Reds", "Oranges", "YlOrBr", "RdPu", "OrRd", "Wistia"]

# Leg 1 of a triplication is solid, leg 2 dashed, leg 3 dotted, etc.,
# cycling if a phase somehow has more legs than this.
LINESTYLES = ["-", "--", ":", "-."]


def _fold_split_indices(x, min_points=6):
    """Split a 1D array into (start, end) index pairs at its local
    extrema, i.e. wherever it stops increasing and starts decreasing
    or vice versa.

    A triplicated phase's distance-vs-time branches are usually
    returned by TauP as one continuous, non-monotonic array (it only
    breaks into separate arrays at true shadow-zone gaps), so this is
    needed to recover the individual legs for separate coloring.
    TauP's internal tau branches (split at every model discontinuity,
    not just at ray turning points) also introduce short spurious
    "legs" a few samples long right at branch boundaries; these are
    merged into the preceding leg rather than kept as their own.

    :param x: Values to split, e.g. a branch's epicentral distances.
    :type x: numpy.ndarray
    :param min_points: Minimum sample count for a run to count as its
        own leg; shorter runs are merged into the previous leg.
    :type min_points: int
    :returns: Inclusive ``(start, end)`` index pairs, one per leg.
    :rtype: list[tuple[int, int]]
    """
    n = len(x)
    if n < 3:
        return [(0, n - 1)] if n else []
    sign = np.sign(np.diff(x))
    for i in range(1, len(sign)):
        if sign[i] == 0:
            sign[i] = sign[i - 1]
    turns = np.where(np.diff(sign) != 0)[0] + 1
    bounds = [0, *turns, n - 1]
    legs = [(bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1)]

    merged = []
    for start, end in legs:
        if merged and (end - start + 1) < min_points:
            merged[-1] = (merged[-1][0], end)
        else:
            merged.append((start, end))
    return merged


def _clip_to_window(dist_km, time_s, max_km):
    """Restrict a monotonic (dist_km, time_s) leg to dist_km <= max_km,
    linearly interpolating one extra point exactly at max_km if the
    leg straddles the window edge.

    Without this, a leg that only has 2 samples (e.g. a head wave like
    Pn/Sn, whose branch is just its critical-distance onset and its
    diffraction cutoff) disappears entirely once the window shrinks to
    where only one of those two points survives a plain boolean mask
    -- matplotlib can't draw a line through a single point. Since a
    head wave's time-distance relation is linear, clipping it exactly
    at the window edge is not an approximation.

    :param dist_km: Monotonic distances for one leg, in km.
    :type dist_km: numpy.ndarray
    :param time_s: Travel times paired with ``dist_km``, in seconds.
    :type time_s: numpy.ndarray
    :param max_km: Window edge, in km.
    :type max_km: float
    :returns: ``(dist_km, time_s)`` clipped to the window.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
    mask = dist_km <= max_km
    if mask.all() or not mask.any():
        return dist_km[mask], time_s[mask]

    idx = np.where(mask)[0]
    lo, hi = idx.min(), idx.max()
    result_d = list(dist_km[idx])
    result_t = list(time_s[idx])
    if hi + 1 < len(dist_km) and not mask[hi + 1]:
        d0, d1 = dist_km[hi], dist_km[hi + 1]
        t0, t1 = time_s[hi], time_s[hi + 1]
        frac = (max_km - d0) / (d1 - d0)
        result_d.append(max_km)
        result_t.append(t0 + frac * (t1 - t0))
    if lo - 1 >= 0 and not mask[lo - 1]:
        d0, d1 = dist_km[lo - 1], dist_km[lo]
        t0, t1 = time_s[lo - 1], time_s[lo]
        frac = (max_km - d0) / (d1 - d0)
        result_d.insert(0, max_km)
        result_t.insert(0, t0 + frac * (t1 - t0))
    return np.array(result_d), np.array(result_t)


def compute_phase_legs(source_depth_km, phase_list=DEFAULT_PHASES,
                        model=DEFAULT_MODEL, max_km=500.0, plot_all=True):
    """Split each requested phase into its individual triplication legs
    and assign each leg its own color.

    This is the shared computation behind :func:`plot_travel_times_km`
    and :func:`plot_ray_paths_km`, so that a given leg (e.g. "Pg"
    turning in the upper vs. lower crust) is drawn in the same color
    on both the travel-time curve and the ray-path fan.

    :param source_depth_km: Source depth, in kilometers.
    :type source_depth_km: float
    :param phase_list: Phase names to include.
    :type phase_list: list[str]
    :param model: Path to an ObsPy TauP ``.npz`` model (or a built-in
        model name, or an existing :class:`~obspy.taup.tau.TauPyModel`).
    :type model: str or obspy.taup.tau.TauPyModel
    :param max_km: Maximum epicentral distance to consider, in km.
    :type max_km: float
    :param plot_all: Also consider the branch mirrored past the
        antipode.
    :type plot_all: bool
    :returns: One dict per leg with keys ``phase``, ``label``,
        ``color``, ``linestyle``, ``seismic_phase`` (the underlying
        :class:`~obspy.taup.seismic_phase.SeismicPhase`), ``lo``/
        ``hi`` (inclusive absolute indices into
        ``seismic_phase.ray_param``/``.dist``/``.time`` for this leg),
        ``dist_km``, ``time_s``, ``mask``, ``mirror_km``, and
        ``mirror_mask``.
    :rtype: list[dict]
    """
    if not isinstance(model, TauPyModel):
        model = TauPyModel(model)

    radius_km = model.model.radius_of_planet
    circumference_km = 2 * np.pi * radius_km
    depth_corrected_model = model.model.depth_correct(source_depth_km)
    phase_names = sorted(parse_phase_list(phase_list))

    cool_i = 0
    warm_i = 0
    phase_cmaps = {}
    legs = []
    for phase in phase_names:
        ph = SeismicPhase(phase, depth_corrected_model)

        # Collect only the legs that actually fall inside the window.
        phase_branches = []
        for s in ph._shadow_zone_splits():
            dist_km_full = ph.dist[s] * radius_km
            time_s_full = ph.time[s]
            for a, b in _fold_split_indices(dist_km_full):
                dist_km = dist_km_full[a:b + 1]
                time_s = time_s_full[a:b + 1]
                mask = dist_km <= max_km
                mirror_km = circumference_km - dist_km
                mirror_mask = mirror_km <= max_km
                if mask.any() or (plot_all and mirror_mask.any()):
                    phase_branches.append((s.start + a, s.start + b,
                                            dist_km, time_s, mask,
                                            mirror_km, mirror_mask))

        n_branches = len(phase_branches)

        # Assign this phase its own colormap: warm for S-family names,
        # cool for everything else (P-family). Different phases within
        # the same family get different colormaps (e.g. Pg -> Blues,
        # Pn -> Greens) so they stay distinguishable from each other,
        # while a phase's own triplication legs are different shades
        # within that one colormap.
        if phase not in phase_cmaps:
            if phase.startswith("S"):
                cmap_name = WARM_CMAPS[warm_i % len(WARM_CMAPS)]
                warm_i += 1
            else:
                cmap_name = COOL_CMAPS[cool_i % len(COOL_CMAPS)]
                cool_i += 1
            phase_cmaps[phase] = plt.get_cmap(cmap_name)
        cmap = phase_cmaps[phase]

        if n_branches <= 1:
            shades = [cmap(0.65)]
        else:
            shades = [cmap(x) for x in np.linspace(0.35, 0.85, n_branches)]

        for j, (lo, hi, dist_km, time_s, mask, mirror_km, mirror_mask) \
                in enumerate(phase_branches):
            color = shades[j]
            linestyle = LINESTYLES[j % len(LINESTYLES)]
            label = phase if n_branches == 1 else f"{phase} (leg {j + 1})"
            legs.append(dict(
                phase=phase, label=label, color=color,
                linestyle=linestyle, seismic_phase=ph,
                lo=lo, hi=hi, dist_km=dist_km, time_s=time_s, mask=mask,
                mirror_km=mirror_km, mirror_mask=mirror_mask,
            ))
    return legs


def plot_travel_times_km(source_depth_km, phase_list=DEFAULT_PHASES,
                          model=DEFAULT_MODEL, max_km=500.0, plot_all=True,
                          legend=True, ax=None, show=True):
    """Plot TauP travel-time curves with distance in kilometers.

    :param source_depth_km: Source depth, in kilometers.
    :type source_depth_km: float
    :param phase_list: Phase names to plot.
    :type phase_list: list[str]
    :param model: Path to an ObsPy TauP ``.npz`` model (or a built-in
        model name, or an existing :class:`~obspy.taup.tau.TauPyModel`).
    :type model: str or obspy.taup.tau.TauPyModel
    :param max_km: Maximum epicentral distance to plot, in kilometers.
    :type max_km: float
    :param plot_all: Also plot the branch mirrored past the antipode.
    :type plot_all: bool
    :param legend: Whether to draw a legend.
    :type legend: bool
    :param ax: Existing axes to plot into.
    :type ax: matplotlib.axes.Axes
    :param show: Whether to call ``plt.show()``.
    :type show: bool
    :returns: The axes used for plotting.
    :rtype: matplotlib.axes.Axes
    """
    legs = compute_phase_legs(source_depth_km, phase_list=phase_list,
                              model=model, max_km=max_km,
                              plot_all=plot_all)

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 12), dpi=100)

    for leg in legs:
        mask = leg["mask"]
        mirror_mask = leg["mirror_mask"]
        plotted = False
        if mask.any():
            d, t = _clip_to_window(leg["dist_km"], leg["time_s"], max_km)
            ax.plot(d, t, label=leg["label"], color=leg["color"],
                    linestyle=leg["linestyle"])
            plotted = True
        if plot_all and mirror_mask.any():
            d, t = _clip_to_window(leg["mirror_km"], leg["time_s"], max_km)
            ax.plot(d, t, label=None if plotted else leg["label"],
                    color=leg["color"], linestyle=leg["linestyle"])

    if legend:
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc="best",
                  numpoints=1)

    ax.grid(True)

    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Time (s)")
    ax.set_xlim(20, 100)
    #ax.set_ylim(bottom=0.0)
    ax.set_ylim(5, 30)

    if show:
        plt.show()
    return ax


def plot_ray_paths_km(source_depth_km, phase_list=DEFAULT_PHASES,
                       model=DEFAULT_MODEL, max_km=500.0,
                       rays_per_leg=6, legend=True, ax=None, show=True):
    """Plot TauP ray paths in a distance-vs-depth cross section, colored
    to match :func:`plot_travel_times_km`'s triplication legs.

    Similar to ``Arrivals.plot_rays(plot_type="cartesian")``, except
    the x-axis is kilometers rather than degrees, and each
    triplication leg (as split by :func:`compute_phase_legs`) is drawn
    in its own color instead of one color per phase name.

    :param source_depth_km: Source depth, in kilometers.
    :type source_depth_km: float
    :param phase_list: Phase names to plot.
    :type phase_list: list[str]
    :param model: Path to an ObsPy TauP ``.npz`` model (or a built-in
        model name, or an existing :class:`~obspy.taup.tau.TauPyModel`).
    :type model: str or obspy.taup.tau.TauPyModel
    :param max_km: Maximum epicentral distance to plot, in kilometers.
    :type max_km: float
    :param rays_per_leg: Number of individual ray paths to draw per
        leg, evenly spaced across its visible distance range.
    :type rays_per_leg: int
    :param legend: Whether to draw a legend.
    :type legend: bool
    :param ax: Existing axes to plot into.
    :type ax: matplotlib.axes.Axes
    :param show: Whether to call ``plt.show()``.
    :type show: bool
    :returns: The axes used for plotting.
    :rtype: matplotlib.axes.Axes
    """
    if not isinstance(model, TauPyModel):
        model = TauPyModel(model)
    radius_km = model.model.radius_of_planet

    legs = compute_phase_legs(source_depth_km, phase_list=phase_list,
                              model=model, max_km=max_km, plot_all=False)

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6), dpi=100)

    for leg in legs:
        ph = leg["seismic_phase"]
        visible_dist_km = leg["dist_km"][leg["mask"]]
        if len(visible_dist_km) == 0:
            continue
        n_rays = min(rays_per_leg, len(visible_dist_km))
        # Evenly spaced target distances across this leg's visible
        # range; each is handed to calc_path, which (unlike shoot_ray)
        # also handles head waves like Pn/Sn.
        targets_km = np.linspace(visible_dist_km.min(),
                                  visible_dist_km.max(), n_rays)

        plotted = False
        for dist_km in targets_km:
            degrees = np.degrees(dist_km / radius_km)
            try:
                arrivals = ph.calc_path(degrees)
            except Exception:
                continue
            # A single distance can have several crossings during a
            # triplication; keep only the one belonging to this leg.
            matches = [a for a in arrivals
                       if leg["lo"] <= a.ray_param_index <= leg["hi"]]
            if not matches:
                continue
            arrival = matches[0]
            path_dist_km = arrival.path["dist"] * radius_km
            path_depth_km = arrival.path["depth"]
            plot_kwargs = dict(label=None if plotted else leg["label"],
                                color=leg["color"],
                                linestyle=leg["linestyle"])
            ax.plot(path_dist_km, path_depth_km, **plot_kwargs)

            # if leg["phase"].startswith("S"):
            #     # Same trick as Arrivals.plot_rays'
            #     # indicate_wave_type=True: matplotlib's "sketch" path
            #     # filter gives S-wave legs a wiggly rather than a
            #     # straight line, so they're distinguishable from P at
            #     # a glance even before checking the legend.
            #     with plt.rc_context({"path.sketch": (2, 10, 1)}):
            #         ax.plot(path_dist_km, path_depth_km, **plot_kwargs)
            # else:
            #     ax.plot(path_dist_km, path_depth_km, **plot_kwargs)
            # plotted = True

    # Reference lines at the model's velocity discontinuities.
    discons = model.model.s_mod.v_mod.get_discontinuity_depths()
    xlim = (0, max_km)
    for depth in discons:
        if depth <= max_km * 2:  # skip absurdly deep ones for a tight plot
            ax.axhline(depth, color="0.75", lw=0.75, zorder=-1)

    ax.plot([0], [source_depth_km], marker="*", color="#FEF215",
            markersize=16, zorder=10, markeredgewidth=1.2,
            markeredgecolor="0.3", clip_on=False)

    if legend:
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc="best",
                  numpoints=1)

    plt.axvline(60)
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Depth (km)")
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 30)
    ax.invert_yaxis()

    if show:
        plt.show()
    return ax


if __name__ == "__main__":
    kwargs = {"max_km": 100,
              "source_depth_km": 0.25
              }
    plot_travel_times_km(show=False, **kwargs)
    plot_ray_paths_km(**kwargs)
