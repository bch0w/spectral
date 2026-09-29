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
# S-family phase its own warm one, sampled once at a fixed point so
# the whole phase (including any triplication branches it folds
# through) is a single consistent color.
COOL_CMAPS = ["Blues", "Greens", "Purples", "BuGn", "GnBu", "PuBu"]
WARM_CMAPS = ["Reds", "OrRd", "YlOrBr", "RdPu", "OrRd", "Wistia"]


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
    """Gather each requested phase's branches and assign each phase its
    own color.

    This is the shared computation behind :func:`plot_travel_times_km`
    and :func:`plot_ray_paths_km`, so that a given phase (e.g. "Sg",
    including whatever triplication branches it folds through) is
    drawn in the same color on both the travel-time curve and the
    ray-path fan. A phase is only split into multiple entries here
    when TauP itself reports a genuine shadow-zone gap
    (``_shadow_zone_splits()``); it is not further split at
    triplication fold points, so a phase retains its own single color
    and linestyle throughout.

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
    :returns: One dict per phase (or per shadow-zone segment, for a
        phase with a genuine gap) with keys ``phase``, ``label``,
        ``color``, ``linestyle``, ``seismic_phase`` (the underlying
        :class:`~obspy.taup.seismic_phase.SeismicPhase`), ``lo``/
        ``hi`` (inclusive absolute indices into
        ``seismic_phase.ray_param``/``.dist``/``.time`` for this
        segment), ``dist_km``, ``time_s``, ``mask``, ``mirror_km``,
        and ``mirror_mask``.
    :rtype: list[dict]
    """
    if not isinstance(model, TauPyModel):
        model = TauPyModel(model)

    radius_km = model.model.radius_of_planet
    circumference_km = 2 * np.pi * radius_km
    depth_corrected_model = model.model.depth_correct(source_depth_km)
    phase_names = sorted(parse_phase_list(phase_list))

    phase_colors = {}
    legs = []
    for i, phase in enumerate(phase_names):
        ph = SeismicPhase(phase, depth_corrected_model)

        # Hardcode look
        if phase == "Pg":
            color = "tab:blue"
            ls = "-"
        elif phase == "Pn":
            color = "tab:green"
            ls="--"
        elif phase == "Sg":
            color = "tab:orange"
            ls="-"
        elif phase == "Sn":
            color = "tab:red"
            ls="--"

        # Collect only the segments that actually fall inside the
        # window.
        for s in ph._shadow_zone_splits():
            dist_km = ph.dist[s] * radius_km
            time_s = ph.time[s]
            mask = dist_km <= max_km
            mirror_km = circumference_km - dist_km
            mirror_mask = mirror_km <= max_km
            if mask.any() or (plot_all and mirror_mask.any()):
                legs.append(dict(
                    phase=phase, label=phase, linestyle=ls, color=color,
                    seismic_phase=ph, lo=s.start, hi=s.start + len(dist_km) - 1,
                    dist_km=dist_km, time_s=time_s, mask=mask,
                    mirror_km=mirror_km, mirror_mask=mirror_mask,
                ))
    return legs


def plot_travel_times_km(source_depth_km, phase_list=DEFAULT_PHASES,
                          model=DEFAULT_MODEL, max_km=500.0, plot_all=True,
                          legend=True, ax=None, show=True, save=False,
                          dpi=100):
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
        _, ax = plt.subplots(figsize=(6, 8), dpi=dpi)

    for leg in legs:
        mask = leg["mask"]
        mirror_mask = leg["mirror_mask"]
        plotted = False
        if mask.any():
            d, t = _clip_to_window(leg["dist_km"], leg["time_s"], max_km)
            ax.plot(d, t, label=leg["label"], color=leg["color"],
                    linestyle=leg["linestyle"], lw=1.75)
            plotted = True
        if plot_all and mirror_mask.any():
            d, t = _clip_to_window(leg["mirror_km"], leg["time_s"], max_km)
            ax.plot(d, t, label=None if plotted else leg["label"],
                    color=leg["color"], linestyle=leg["linestyle"], lw=1.75)

    if legend:
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc="best",
                  numpoints=1, fontsize=15, frameon=False)

    ax.grid(True)

    ax.set_xlabel("Distance (km)", fontsize=16)
    ax.set_ylabel("Time (s)", fontsize=16)
    ax.set_title("ak135f travel time curve", fontsize=16)
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(1.25)
    ax.set_xlim(20, 100)
    #ax.set_ylim(bottom=0.0)
    ax.set_ylim(5, 25)

    plt.axvline(60, c="k", zorder=3, lw=2, ls="--")

    # Without this, the y-axis label (and the tick labels next to it)
    # can get clipped at the left edge of the figure, since the axes
    # box doesn't otherwise leave room for a fontsize=16 label.
    ax.figure.tight_layout()

    if save:
        plt.savefig(save, bbox_inches="tight")
    if show:
        plt.show()

    return ax


def plot_ray_paths_km(source_depth_km, phase_list=DEFAULT_PHASES,
                       model=DEFAULT_MODEL, max_km=500.0,
                       wiggly_s=True, legend=True, ax=None, show=True,
                       save=False, dpi=100):
    """Plot TauP ray paths in a distance-vs-depth cross section, colored
    to match :func:`plot_travel_times_km`'s phase colors.

    Similar to ``Arrivals.plot_rays(plot_type="cartesian")``, except
    the x-axis is kilometers rather than degrees, and colors are
    assigned per :func:`compute_phase_legs` (cool colors for P-family
    phases, warm for S-family).

    Each phase contributes exactly one ray: the one that arrives at
    ``max_km``. A phase whose reachable distance range doesn't extend
    out to ``max_km`` is skipped rather than drawing a shorter ray, so
    every line in the plot ends exactly at the right edge instead of
    stopping short of it or, if picked from some other distance,
    running past it.

    :param source_depth_km: Source depth, in kilometers.
    :type source_depth_km: float
    :param phase_list: Phase names to plot.
    :type phase_list: list[str]
    :param model: Path to an ObsPy TauP ``.npz`` model (or a built-in
        model name, or an existing :class:`~obspy.taup.tau.TauPyModel`).
    :type model: str or obspy.taup.tau.TauPyModel
    :param max_km: Epicentral distance each ray path should arrive at,
        in kilometers; also the plot's x-axis limit.
    :type max_km: float
    :param wiggly_s: If ``True``, draw legs whose phase name starts
        with "S" as a wiggly line via matplotlib's "sketch" path
        filter (the same trick ``Arrivals.plot_rays`` uses with
        ``indicate_wave_type=True``), so S-wave legs are visually
        distinguishable from P-wave legs even without checking colors
        or the legend. Set to ``False`` for plain straight/dashed
        lines throughout.
    :type wiggly_s: bool
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
        _, ax = plt.subplots(figsize=(6, 8), dpi=dpi)

    degrees = np.degrees(max_km / radius_km)
    for leg in legs:
        ph = leg["seismic_phase"]
        # Skip legs that don't actually reach out to max_km -- there's
        # no arrival to draw for them at this distance.
        if not (leg["dist_km"].min() <= max_km <= leg["dist_km"].max()):
            continue
        try:
            arrivals = ph.calc_path(degrees)
        except Exception:
            continue
        # A single distance can have several crossings during a
        # triplication; keep only ones belonging to this phase segment
        # (relevant if the phase has more than one shadow-zone
        # segment) and just take the first -- with legs no longer
        # split at fold points, any of them is an equally valid single
        # representative ray for this phase.
        matches = [a for a in arrivals
                   if leg["lo"] <= a.ray_param_index <= leg["hi"]]
        if not matches:
            continue
        arrival = matches[0]
        for arrival in arrivals:
            path_dist_km = arrival.path["dist"] * radius_km
            path_depth_km = arrival.path["depth"]
            plot_kwargs = dict(label=leg["label"], color=leg["color"])
            if wiggly_s and leg["phase"].startswith("S"):
                with plt.rc_context({"path.sketch": (6, 20, 1)}):
                    ax.plot(path_dist_km, path_depth_km, lw=1, **plot_kwargs)
            else:
                ax.plot(path_dist_km, path_depth_km, linestyle=leg["linestyle"],
                        zorder=10, **plot_kwargs)

    # Reference lines at the model's velocity discontinuities.
    discons = model.model.s_mod.v_mod.get_discontinuity_depths()
    for depth in discons:
        if depth <= max_km * 2:  # skip absurdly deep ones for a tight plot
            ax.axhline(depth, color="0.75", lw=2, zorder=-1)

    ax.plot([0], [source_depth_km], marker="*", color="#FEF215",
            markersize=16, zorder=10, markeredgewidth=1.2,
            markeredgecolor="k", clip_on=False)

    if legend:
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc="lower right",
                  numpoints=1, frameon=False, fontsize=15)

    ax.set_xlabel("Distance (km)", fontsize=16)
    ax.set_ylabel("Depth (km)", fontsize=16)
    ax.set_title("ak135f ray paths", fontsize=16)
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(1.25)
    ax.set_xlim(-2, max_km+1)
    ax.set_ylim(-2, 21)

    # Annotate each layer's Vp/Vs, right at the depth where it ends.
    # A layer that extends past the bottom of the plot (the "mantle"
    # layer starting at 18 km here, whose true base is at 80 km) gets
    # its annotation placed just above the lower y limit instead, and
    # its velocities are evaluated at that same clipped depth rather
    # than at its real, off-screen bottom.
    v_mod = model.model.s_mod.v_mod
    y_top, y_bottom = sorted(ax.get_ylim())
    margin = (y_bottom - y_top) * 0.02
    x_text = max_km * 0.02
    for top, bot in zip(discons[:-1], discons[1:]):
        if top >= y_bottom:
            break
        y_text = bot if bot <= y_bottom else y_bottom - margin
        vp = v_mod.evaluate_above(y_text, "p").item()
        vs = v_mod.evaluate_above(y_text, "s").item()
        ax.text(x_text, y_text, f"Vp={vp:.2f}, Vs={vs:.2f}", fontsize=12,
                color="k", va="bottom", ha="left", zorder=10,
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none",
                                          pad=1))
        if bot > y_bottom:
            break

    ax.invert_yaxis()

    # Same fix as plot_travel_times_km: without this, the y-axis label
    # (and its tick labels) can get clipped at the left edge of the
    # figure.
    ax.figure.tight_layout()

    if save:
        plt.savefig(save, bbox_inches="tight")
    if show:
        plt.show()

    return ax


if __name__ == "__main__":
    kwargs = {"dpi": 250,
              "source_depth_km": 0.25,
              "phase_list": ["Pg", "Pn", "Sg", "Sn"], 
              }
    plot_travel_times_km(max_km=120, show=False, save="ak135f_tt.png", **kwargs)
    plot_ray_paths_km(max_km=60, save="ak135f_rp.png", **kwargs)
