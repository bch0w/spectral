"""
Distance Record Section (DistRecSec)

For a single, fixed source mechanism, query Instaseis/AxiSEM synthetics at a
range of source-receiver distances (fixed backazimuth, i.e. a great-circle
profile) and plot them as a true record section: distance on the y-axis,
time on the x-axis. Theoretical TauP arrival-time curves are overlaid for
user-defined phases so that phase moveout, crossover distances (e.g.
Pn/Pg, Sn/Sg) and phase conversions can be tracked visually as a function of
distance from the source.

Unlike `momtenmeas.py`'s record sections (which offset waveforms by a fixed,
arbitrary spacing and mark single-distance TauP picks as vertical lines),
here waveforms are offset by their *actual* distance and TauP arrivals are
plotted as continuous distance-vs-time curves, since the whole point is to
watch how a phase's timing (and existence) changes across the profile.
"""
import hashlib
import os
import sys
import time
import instaseis
import matplotlib.pyplot as plt
import numpy as np

from concurrent.futures import ProcessPoolExecutor, as_completed
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from obspy import read, read_events, Stream
from obspy.geodetics import kilometers2degrees, gps2dist_azimuth
from obspy.taup import TauPyModel
from scipy.ndimage import maximum_filter1d, uniform_filter1d

# Import plotting helpers from the internal spectral Class
sys.path.insert(0, "/Users/prof/Repos/spectral")
from plotools.prettyplot import cmaphex, set_plot_aesthetic


# Default TauP phase list to look for evolution/conversion of body waves
DEFAULT_PHASES = ["Pn", "Pg", "Sn", "Sg"]


class DistRecSec:
    """
    For a fixed source mechanism, query Instaseis/AxiSEM synthetics at a list
    of source-receiver distances and cache them to disk so that they can be
    reloaded (and re-plotted) without re-querying Instaseis each time.
    """
    def __init__(self, dist_km_list, baz, src_depth_km, tmin=None, tmax=None,
                 corners=4, components="Z", kind="velocity",
                 syngine="iasp91_2s", taup_model="iasp91",
                 taper_vel_kmps=None, taper_vel_width_s=5.,
                 fig_path="FIGURES", wav_path="SAC"):
        """
        Set up calling structure

        Parameters
        ----------
        dist_km_list : list of float
            Source-receiver distances in km to gather synthetics for
        baz : float
            Backazimuth from source to receivers in degrees. Fixed for all
            distances so that all stations sit on the same great-circle
            profile
        src_depth_km : float
            Source depth in kilometers
        tmin : float, optional
            Minimum filter period in seconds. Waveforms are filtered and
            saved filtered so that plotting always reflects what was
            processed (see `momtenmeas.py`)
        tmax : float, optional
            Maximum filter period in seconds
        corners : int, optional
            Number of filter corners, default 4
        components : str, optional
            Seismogram components to extract, default is "Z"
        kind : str, optional
            Waveform type ('velocity', 'displacement', etc.), default is
            "velocity"
        syngine : str, optional
            Syngine model to use for synthetics, default is "iasp91_2s"
        taup_model : str, optional
            TauP velocity model matching `syngine`, used to draw arrival
            curves, default is "iasp91"
        taper_vel_kmps : float, optional
            If given, mute (taper to 0) all signal arriving after this group
            velocity [km/s] for each station, based on its source-receiver
            distance -- e.g. 2.5 removes everything slower than 2.5 km/s
            (typically surface waves/coda). See `taper_after_velocity()`.
            Applied during `process()`, alongside the bandpass filter.
            Default None (no taper)
        taper_vel_width_s : float, optional
            Taper duration in seconds used to smoothly ramp to zero past the
            `taper_vel_kmps` cutoff (avoids ringing from an abrupt mute),
            default 5.0
        fig_path : str, optional
            Directory path to save figures, default is "FIGURES"
        wav_path : str, optional
            Directory path to save waveform data, default is "SAC"
        """
        self.dist_km_list = sorted(dist_km_list)
        self.baz = baz
        self.src_depth_km = src_depth_km
        self.tmin = tmin
        self.tmax = tmax
        self.corners = corners
        self.components = components
        self.kind = kind
        self.taper_vel_kmps = taper_vel_kmps
        self.taper_vel_width_s = taper_vel_width_s

        self.syngine = syngine
        print(f"opening instaseis db 'syngine://{self.syngine}'")
        self.db = instaseis.open_db(f"syngine://{self.syngine}")
        self.taup_model = taup_model

        # !!! BCBC (see `momtenmeas.py`): Syngine strips the water/mud layer
        # so the top of the model is solid Earth; source is nudged just below
        # the surface to avoid sitting exactly on it
        self.model_surface_m = 0

        self.fig_path = fig_path
        self.wav_path = wav_path
        for path_ in [self.fig_path, self.wav_path]:
            os.makedirs(path_, exist_ok=True)

        # Populated by `get_src()` / `run()`
        self.src = None
        self.src_tag = None
        self.save_tag = None

        print(f"DistRecSec set up for {len(self.dist_km_list)} distances "
              f"[{min(self.dist_km_list):.1f}, {max(self.dist_km_list):.1f}]km; "
              f"baz={self.baz}deg; z={self.src_depth_km}km")

    def get_src(self, path_cmtsolution, mt=None, latitude=0, longitude=0):
        """
        Build the single, fixed Instaseis Source used for the whole record
        section. Provide exactly one of:

        - `path_cmtsolution`: path to a CMTSOLUTION/QuakeML file, moment
          tensor read via ObsPy `read_events` (as in `momtenmeas.py`)
        - `mt`: dict with keys m_rr, m_tt, m_pp, m_rt, m_rp, m_tp [Nm]
        """
        depth_in_m = self.src_depth_km * 1e3 + self.model_surface_m

        print(f"reading source mechanism from '{path_cmtsolution}'")
        event = read_events(path_cmtsolution)[0]
        tensor = event.focal_mechanisms[0].moment_tensor.tensor
        mt = {"m_rr": tensor["m_rr"], "m_tt": tensor["m_tt"],
              "m_pp": tensor["m_pp"], "m_rt": tensor["m_rt"],
              "m_rp": tensor["m_rp"], "m_tp": tensor["m_tp"]}

        print(f"building source from moment tensor: {mt}")
        self.src = instaseis.Source(
            latitude=latitude, longitude=longitude,
            depth_in_m=depth_in_m, **mt
            )

        print(f"\tsource: lat={self.src.latitude:.2f} "
              f"lon={self.src.longitude:.2f} "
              f"depth={self.src.depth_in_m * 1e-3:.2f}km")

        # Deterministic tag identifying this exact mechanism, used below in
        # `run()` to key the on-disk SAC cache. Derived from `self.src`'s
        # actually-resolved tensor components (not from which input branch
        # was used above) so it's correct regardless of whether the caller
        # went through `path_cmtsolution`, `mt`, or `strike`/`dip`/`rake` --
        # and so a caller that skips `main()` and drives `DistRecSec`
        # directly still can't collide two different mechanisms onto the
        # same cache files
        if path_cmtsolution is not None:
            self.src_tag = os.path.basename(str(path_cmtsolution))
        else:
            vals = "_".join(f"{getattr(self.src, k):.6g}" for k in
                            ["m_rr", "m_tt", "m_pp", "m_rt", "m_rp", "m_tp"])
            self.src_tag = "mt" + hashlib.md5(vals.encode()).hexdigest()[:8]
        print(f"\tsrc_tag='{self.src_tag}'")

        return self.src

    def get_rcv(self, dist_km, network="XX", station=None):
        """
        Define a receiver at `dist_km` from the source along the fixed
        backazimuth `self.baz`. Rudimentary Cartesian shooting, same
        approach as `momtenmeas.py`, fine for regional/teleseismic profiles
        """
        # Instaseis enforces station codes of <= 5 chars, so round to int km
        station = station or f"{int(round(dist_km)):0>5d}"
        dist_deg = kilometers2degrees(dist_km)

        lon = dist_deg * np.cos(np.deg2rad(self.baz)) + self.src.longitude
        lat = dist_deg * np.sin(np.deg2rad(self.baz)) + self.src.latitude

        return instaseis.Receiver(latitude=lat, longitude=lon,
                                  network=network, station=station)

    def get_synthetics(self, dist_km):
        """
        Query Instaseis for a single `dist_km`, attach SAC headers and save
        the raw (unfiltered) waveform to disk so repeat calls don't requery
        """
        print(f"\tquerying instaseis: d={dist_km:.1f}km "
              f"comp={self.components} kind={self.kind}")
        rcv = self.get_rcv(dist_km)
        st = self.db.get_seismograms(source=self.src, receiver=rcv,
                                     components=self.components,
                                     kind=self.kind, dt=0.1,
                                     remove_source_shift=True)

        dist_m, az, baz = gps2dist_azimuth(lat1=self.src.latitude,
                                           lon1=self.src.longitude,
                                           lat2=rcv.latitude, lon2=rcv.longitude)

        sac_header = {
            "iztype": 9,  # Ref time equivalence, IB (9): Begin time
            "b": st[0].stats.starttime,
            "e": st[0].stats.endtime,
            "evla": self.src.latitude,
            "evlo": self.src.longitude,
            "evdp": self.src.depth_in_m * 1e-3,
            "stla": rcv.latitude,
            "stlo": rcv.longitude,
            "stel": 0,
            "kevnm": str(self.save_tag)[:16],
            "nzyear": st[0].stats.starttime.year,
            "nzjday": st[0].stats.starttime.julday,
            "nzhour": st[0].stats.starttime.hour,
            "nzmin": st[0].stats.starttime.minute,
            "nzsec": st[0].stats.starttime.second,
            "nzmsec": int(f"{st[0].stats.starttime.microsecond:0>6}"[:3]),
            "dist": dist_m * 1e-3,
            "az": az,
            "baz": baz,
            "gcarc": kilometers2degrees(dist_m * 1e-3),
            "lpspol": 0,
            "lcalda": 1,
        }

        for tr in st:
            tr.stats.sac = sac_header
            fid = (f"{self.wav_path}/{self.save_tag}_d{dist_km:g}_"
                   f"{tr.stats.component}.SAC")
            tr.write(fid, format="SAC")
        print(f"\t\tsaved {len(st)} trace(s) -> "
              f"{self.wav_path}/{self.save_tag}_d{dist_km:g}_*.SAC")

        return st

    def load_synthetics(self, dist_km):
        """
        Load previously saved SAC files for `dist_km` from disk, raises if
        any component is missing so caller can fall back to `get_synthetics`
        """
        st = Stream()
        for comp in self.components:
            fid = (f"{self.wav_path}/{self.save_tag}_d{dist_km:g}_"
                   f"{comp}.SAC")
            if os.path.exists(fid):
                st += read(fid)
            else:
                raise FileNotFoundError(fid)
        print(f"\tloaded cached synthetics: d={dist_km:.1f}km "
              f"({self.components})")
        return st

    def process(self, dist_km):
        """
        Get (querying Instaseis if not already saved), filter, and
        group-velocity-taper synthetics for a single `dist_km`. Split out
        from `run()` so it can be mapped over `dist_km_list`, including in
        parallel
        """
        try:
            st = self.load_synthetics(dist_km)
        except FileNotFoundError:
            st = self.get_synthetics(dist_km)

        if self.tmin and self.tmax:
            print(f"\tfiltering d={dist_km:.1f}km "
                  f"[{self.tmin}, {self.tmax}]s corners={self.corners}")
            st.filter("bandpass", freqmin=1 / self.tmax, freqmax=1 / self.tmin,
                      zerophase=False, corners=self.corners)

        if self.taper_vel_kmps:
            t_cut = dist_km / self.taper_vel_kmps
            print(f"\ttapering d={dist_km:.1f}km after "
                  f"{self.taper_vel_kmps}km/s (t={t_cut:.1f}s, "
                  f"width={self.taper_vel_width_s}s)")
            for tr in st:
                tr.data = taper_after_velocity(
                    tr.data.astype(float), tr.times(), dist_km,
                    self.taper_vel_kmps, taper_s=self.taper_vel_width_s
                    )

        return dist_km, st

    def run(self, save_tag="recsec", parallel=True):
        """
        Loop over `self.dist_km_list`, get/load + filter synthetics for each
        distance. Requires `get_src()` to have been called first.

        `save_tag` is a caller-supplied context label (e.g. syngine/depth/
        baz); the actual cache-file tag always appends `self.src_tag` (set
        by `get_src()`) on top of it, so this method can never silently
        collide two different source mechanisms onto the same cached SAC
        files even if the caller reuses the same `save_tag` string

        Returns
        -------
        dict of {dist_km: Stream}
        """
        assert self.src is not None, "call `get_src()` before `run()`"
        self.save_tag = f"{save_tag}_{self.src_tag}"

        n = len(self.dist_km_list)
        print(f"gathering synthetics for {n} distances "
              f"(save_tag='{self.save_tag}', parallel={parallel})")
        t0 = time.time()

        dist_dict = {}
        if parallel:
            with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
                futures = {executor.submit(self.process, dist_km): dist_km
                          for dist_km in self.dist_km_list}
                for i, future in enumerate(as_completed(futures), 1):
                    dist_km, st = future.result()
                    dist_dict[dist_km] = st
                    print(f"[{i}/{n}] done d={dist_km:.1f}km")
        else:
            for i, dist_km in enumerate(self.dist_km_list, 1):
                dist_km, st = self.process(dist_km)
                dist_dict[dist_km] = st
                print(f"[{i}/{n}] done d={dist_km:.1f}km")

        print(f"finished gathering {n} distances in "
              f"{time.time() - t0:.1f}s")

        return dist_dict


def taup_phase_curves(dist_km_list, src_depth_km, tp_phases,
                      taup_model="iasp91", dist_step_km=None, n_dist=300):
    """
    Compute TauP arrival-time curves as a function of distance for each phase
    in `tp_phases`, sampling a dense array of distances (independent of the
    actual station spacing) so curves are smooth and show exactly where a
    phase appears, disappears or is overtaken by another (crossover/
    conversion distances).

    Parameters
    ----------
    dist_km_list : list of float
        Used only to get the [min, max] distance range to sample over
    src_depth_km : float
        Source depth in km
    tp_phases : list of str
        TauP phase names to query, e.g. ["Pn", "Pg", "Sn", "Sg"]
    taup_model : str, optional
        TauP velocity model, default "iasp91"
    dist_step_km : float, optional
        If given, sample the distance range in fixed `dist_step_km`
        increments instead of `n_dist` evenly spaced points
    n_dist : int, optional
        Number of distances to sample if `dist_step_km` not given

    Returns
    -------
    dict of {phase_name: (dist_km_array, time_s_array)}
        Only the first (fastest) arrival of a given phase name is kept at
        each distance so each phase maps to a single, clean curve. Distances
        where the phase does not exist are simply omitted (not NaN-padded),
        which is what creates the visual break at a crossover/cutoff
        distance when plotted with `ax.plot`
    """
    model = TauPyModel(model=taup_model)
    dmin, dmax = min(dist_km_list), max(dist_km_list)
    if dist_step_km:
        dist_arr = np.arange(dmin, dmax + dist_step_km, dist_step_km)
    else:
        dist_arr = np.linspace(dmin, dmax, n_dist)

    print(f"computing TauP arrival curves: {len(tp_phases)} phases "
          f"({', '.join(tp_phases)}) across {len(dist_arr)} distances "
          f"[{dmin:.1f}, {dmax:.1f}]km, model='{taup_model}'")
    t0 = time.time()

    out = {phase: {"dist": [], "time": []} for phase in tp_phases}
    for dist_km in dist_arr:
        dist_deg = kilometers2degrees(dist_km)
        arrivals = model.get_travel_times(source_depth_in_km=src_depth_km,
                                          distance_in_degree=dist_deg,
                                          phase_list=tp_phases)
        seen = set()
        for arrival in arrivals:
            if arrival.name in seen:
                continue
            seen.add(arrival.name)
            out[arrival.name]["dist"].append(dist_km)
            out[arrival.name]["time"].append(arrival.time)

    curves = {phase: (np.array(v["dist"]), np.array(v["time"]))
             for phase, v in out.items()}

    for phase, (d, t) in curves.items():
        if len(d):
            print(f"\t{phase}: present {d.min():.1f}-{d.max():.1f}km "
                  f"({len(d)} pts)")
        else:
            print(f"\t{phase}: not found across this distance range")
    print(f"\tdone in {time.time() - t0:.1f}s")

    return curves


def agc(data, dt, window_s, method="rms", eps_pct=1e-2, return_gain=False):
    """
    Automatic Gain Control (AGC): rescale a trace by a sliding local
    amplitude estimate so that weak/early arrivals become visible at a
    scale comparable to the dominant (usually much larger) arrival, instead
    of being flattened by a single whole-trace normalization. This is the
    standard trick used in seismic record sections for exactly this problem.

    This is a *display-only* transform -- it operates on and returns a
    plain numpy array, never touching the underlying ObsPy Trace/Stream, so
    the actual waveform data is never modified.

    Parameters
    ----------
    data : np.ndarray
        Waveform samples
    dt : float
        Sample spacing in seconds
    window_s : float
        Sliding window length in seconds used to estimate the local gain.
        Shorter windows equalize amplitude more aggressively (every wiggle
        ends up a similar size); longer windows are gentler and preserve
        more of the true relative amplitude between nearby arrivals
    method : str, optional
        'rms' (root-mean-square amplitude in the window, smoother) or 'max'
        (peak absolute amplitude in the window, punchier), default 'rms'
    eps_pct : float, optional
        Stabilizes the gain in low-energy regions (e.g. before the first
        arrival) so near-silence doesn't get amplified into noise.
        Expressed as a fraction of the trace's overall RMS/max, default 1e-2
    return_gain : bool, optional
        If True, also return the per-sample multiplicative gain factor that
        was applied (1 / (envelope + eps)) -- useful for e.g. color-coding a
        plot by how much a given sample was boosted, default False

    Returns
    -------
    np.ndarray
        Gain-corrected trace, same length as `data`
    np.ndarray, optional
        Per-sample gain factor applied, only if `return_gain=True`
    """
    data = data.astype(float)
    win = max(int(round(window_s / dt)), 1)

    if method == "rms":
        # clip: uniform_filter1d can return tiny negative values for
        # near-zero-energy windows due to floating-point roundoff
        mean_sq = np.clip(uniform_filter1d(data ** 2, size=win), 0, None)
        envelope = np.sqrt(mean_sq)
        eps = eps_pct * np.sqrt(np.mean(data ** 2))
    elif method == "max":
        envelope = maximum_filter1d(np.abs(data), size=win)
        eps = eps_pct * np.max(np.abs(data))
    else:
        raise ValueError("`method` must be 'rms' or 'max'")

    gain = 1. / (envelope + eps)
    gained = data * gain

    if return_gain:
        return gained, gain
    return gained


def taper_after_velocity(data, times, dist_km, vel_kmps, taper_s=5.):
    """
    Mute (taper to 0) all signal arriving after a given group velocity, for
    a trace at a known source-receiver distance. Commonly used to remove
    late-arriving energy (e.g. surface waves/coda past the body-wave window)
    so it doesn't contaminate analysis or dominate a record-section display.

    The cutoff time is `t_cut = dist_km / vel_kmps`; a smooth half-cosine
    (raised-cosine) taper ramps the signal to 0 over `[t_cut, t_cut +
    taper_s]` rather than cutting it off abruptly, which would ring.

    Parameters
    ----------
    data : np.ndarray
        Waveform samples
    times : np.ndarray
        Time values [s] for each sample in `data`, relative to the same
        origin time used to compute `dist_km / vel_kmps` (i.e. the trace's
        own `tr.times()`)
    dist_km : float
        Source-receiver distance in km for this trace
    vel_kmps : float
        Group velocity [km/s]; signal arriving after `dist_km / vel_kmps`
        seconds is tapered out
    taper_s : float, optional
        Taper duration in seconds from `t_cut` to fully zeroed, default 5.0

    Returns
    -------
    np.ndarray
        Tapered copy of `data`, same length
    """
    t_cut = dist_km / vel_kmps
    mask = np.ones_like(data, dtype=float)

    if taper_s > 0:
        ramp = (times >= t_cut) & (times < t_cut + taper_s)
        mask[ramp] = 0.5 * (1 + np.cos(np.pi * (times[ramp] - t_cut) / taper_s))
    mask[times >= t_cut + taper_s] = 0.

    return data * mask


def reduce_time_shift_fn(dists, reduce_by, src_depth_km=None,
                         taup_model="iasp91", dist_step_km=None, n_dist=300):
    """
    Build a function that gives the time shift [s] to subtract from a
    trace at a given distance, for a reduced-velocity or phase-aligned
    record section.

    Parameters
    ----------
    dists : list of float
        Distances (km) spanning the record section; only used to set the
        distance range of the reference curve when `reduce_by` is a phase
        name
    reduce_by : float or str
        - float: a constant reduction velocity [km/s]. Classic "reduced
          travel-time" record section -- shift = dist_km / reduce_by. A
          phase moving at exactly this apparent velocity becomes flat
          (horizontal) across the section; other phases are only
          approximately flattened
        - str: a TauP phase name (e.g. 'Pg'). Every trace is shifted by the
          *actual* TauP first-arrival time of that phase at its own
          distance (via `taup_phase_curves()`), so traces are exactly
          aligned on that phase's arrival rather than just approximately
          flattened by a single velocity
    src_depth_km : float
        Source depth in km, required when `reduce_by` is a phase name
    taup_model : str, optional
        TauP velocity model, required when `reduce_by` is a phase name
    dist_step_km / n_dist : see `taup_phase_curves`

    Returns
    -------
    callable
        `shift_fn(dist_km)` -> time shift in seconds (accepts a scalar or
        an array of distances)
    """
    if isinstance(reduce_by, str):
        curves = taup_phase_curves(dists, src_depth_km, [reduce_by],
                                   taup_model=taup_model,
                                   dist_step_km=dist_step_km, n_dist=n_dist)
        ref_dist, ref_time = curves[reduce_by]
        if len(ref_dist) == 0:
            raise ValueError(f"no '{reduce_by}' arrivals found across "
                             f"[{min(dists):.1f}, {max(dists):.1f}]km to "
                             f"align on")
        # Beyond the phase's own distance range, np.interp holds the
        # nearest edge value (flat extrapolation) rather than erroring
        return lambda d: np.interp(d, ref_dist, ref_time)

    return lambda d: np.asarray(d) / reduce_by


def plot_record_section(dist_dict, src_depth_km, tp_phases=None,
                        taup_model="iasp91", component="Z", wf_scale=1.,
                        wf_color="k", wf_lw=0.8, fill=True, reduce_by=None,
                        agc_window_s=None, agc_method="rms", agc_eps_pct=1e-2,
                        gain_cmap="inferno", gain_clip_pct=(1, 99),
                        dist_step_km=None, n_dist=300, tp_cmap="rainbow_r",
                        tp_lw=1.5, xlim=None, ylim=None, title=None,
                        xlabel=None, ylabel="Distance [km]",
                        ytick_step_km=None, ytick_fontsize=8,
                        fig_size=(10, 12), dpi=100, save=None, show=True):
    """
    Plot a distance record section: waveforms plotted at y = distance [km],
    scaled by `wf_scale`. TauP arrival-time curves for `tp_phases` are
    overlaid across the full distance range to track phase moveout/
    conversion with distance.

    By default all traces share a single global normalization (the peak
    absolute amplitude across the whole `dist_dict`), so relative amplitude
    differences between distances (e.g. geometrical spreading/attenuation)
    are preserved in the display rather than washed out by normalizing each
    trace to its own max independently. Because that necessarily means
    far/small traces can be dwarfed by near/large ones -- and even within a
    single trace, a dominant arrival (e.g. Pg) can visually swamp smaller,
    earlier or later phases (e.g. Pn, or converted phases) -- set
    `agc_window_s` to apply Automatic Gain Control: a sliding local-
    amplitude equalization (see `agc()`) that makes weak arrivals visible at
    a comparable scale to the dominant one. AGC is a display-only transform;
    it never modifies the underlying Stream data.

    When AGC is active, each waveform is drawn as a `LineCollection` colored
    sample-by-sample according to how much gain was applied at that point
    (in dB, `gain_cmap` colormap, dark = little/no gain = signal was already
    close to full scale, bright = lots of gain = signal was much smaller
    than the loudest arrival before boosting), with a shared colorbar so
    color is directly comparable across the whole record section -- both
    between distances (far/weak stations read uniformly brighter) and within
    a single trace (a quiet early phase reads brighter than the dominant
    arrival next to it).

    Parameters
    ----------
    dist_dict : dict
        {dist_km: Stream}, output of `DistRecSec.run()`
    src_depth_km : float
        Source depth in km, required to compute TauP curves
    tp_phases : list of str, optional
        TauP phase names to overlay as arrival-time curves. If not given, no
        curves are plotted
    taup_model : str, optional
        TauP velocity model, should match the one used to generate synthetics
    component : str, optional
        Which component to plot from each Stream, default "Z"
    wf_scale : float, optional
        Amplitude scaling for the (normalized) waveforms, in the same units
        as the y-axis (km), controls how much waveforms overlap
    wf_color : str, optional
        Waveform fill color (and line color when `agc_window_s` is not set
        or `gain_cmap` is None)
    wf_lw : float, optional
        Waveform linewidth
    fill : bool, optional
        Fill positive wiggle excursions for a classic record-section look
    reduce_by : float or str, optional
        Shift each trace's time axis by an amount that depends on its
        distance, for a reduced-velocity / phase-aligned record section
        (see `reduce_time_shift_fn()`). This is a display-only shift, both
        the waveforms and the TauP overlay curves are shifted consistently,
        the underlying Stream data is never modified:
        - float: constant reduction velocity [km/s] (classic "reduced
          travel-time" plot -- shift = dist_km / reduce_by). A phase moving
          at exactly this apparent velocity becomes flat/horizontal
        - str: a TauP phase name, e.g. 'Pg' -- every trace is shifted by
          the *actual* TauP arrival time of that phase at its own distance,
          so all traces align exactly on that phase (t=0), which is more
          precise than a constant-velocity reduction since it accounts for
          the phase's true (non-linear) moveout curve
        Default None (absolute time axis, no shift)
    agc_window_s : float, optional
        If given, apply Automatic Gain Control with this sliding window
        length in seconds (see `agc()`) to each trace before plotting, so
        that weak/early arrivals aren't swamped by the dominant arrival.
        If not given (default), traces are shown at true relative amplitude
        (global normalization only, no local gain correction)
    agc_method : str, optional
        'rms' or 'max', passed to `agc()`, default 'rms'
    agc_eps_pct : float, optional
        Stabilization factor passed to `agc()`, default 1e-2
    gain_cmap : str, optional
        Colormap used to color each waveform by the AGC gain applied at
        each sample (only used when `agc_window_s` is set). Set to None to
        keep plain `wf_color` lines even with AGC on. Default 'inferno'
        (dark -> bright), a good match for "little gain" -> "lots of gain"
    gain_clip_pct : (float, float), optional
        Percentile range used to set the color scale limits from the pooled
        gain (in dB) across all traces, default (1, 99). Prevents a handful
        of near-silent, extreme-gain samples (e.g. before the first arrival)
        from washing out the color contrast over the rest of the section
    dist_step_km / n_dist : see `taup_phase_curves`
    tp_cmap : str, optional
        Colormap used to give each TauP phase a distinct color
    tp_lw : float, optional
        Linewidth for TauP arrival curves
    xlim / ylim : list, optional
        Axis limits, ylim in km, xlim in s
    title : str, optional
    xlabel / ylabel : str, optional
    ytick_step_km : float, optional
        Label the y-axis at fixed `ytick_step_km` increments instead of at
        each individual station distance. Useful to de-clutter dense
        record sections with many closely-spaced distances
    ytick_fontsize : float, optional
        Font size for the distance tick labels
    fig_size : tuple, optional
    dpi : float, optional
    save : str, optional
        Filename to save figure to
    show : bool, optional

    Returns
    -------
    (Figure, Axes)
    """
    dists = sorted(dist_dict.keys())
    traces = {d: dist_dict[d].select(component=component)[0] for d in dists}

    # Global (not per-trace) normalization so relative amplitude between
    # distances stays physically meaningful
    global_max = max(np.max(np.abs(tr.data)) for tr in traces.values()) or 1.

    # Reduced-velocity / phase-aligned time axis: a per-distance shift
    # subtracted from both the waveforms and the TauP overlay curves below
    shift_fn = None
    if reduce_by is not None:
        shift_fn = reduce_time_shift_fn(dists, reduce_by,
                                        src_depth_km=src_depth_km,
                                        taup_model=taup_model,
                                        dist_step_km=dist_step_km,
                                        n_dist=n_dist)
        shifts = shift_fn(np.array(dists))
        print(f"\treducing time axis by '{reduce_by}': shifts range "
              f"[{shifts.min():.1f}, {shifts.max():.1f}]s across "
              f"[{min(dists):.1f}, {max(dists):.1f}]km")

    if xlabel is None:
        if reduce_by is None:
            xlabel = "Time [s]"
        elif isinstance(reduce_by, str):
            xlabel = f"Time - {reduce_by} arrival [s]"
        else:
            xlabel = f"Reduced Time [t - x/{reduce_by:g}km/s] (s)"

    print(f"plotting record section: {len(dists)} distances, "
          f"comp={component}, wf_scale={wf_scale}, global_max={global_max:.3E}, "
          f"agc_window_s={agc_window_s}, reduce_by={reduce_by}")
    t0 = time.time()

    f, ax = plt.subplots(figsize=fig_size, dpi=dpi)

    # First pass: compute the display waveform (and AGC gain, if requested)
    # for every distance, so a single, shared color scale can be built
    # before any plotting happens
    processed = {}
    gain_db_all = []
    for dist_km in dists:
        tr = traces[dist_km]
        raw = tr.data.astype(float)

        if agc_window_s:
            # AGC equalizes local amplitude to ~unit scale itself, so scale
            # directly by `wf_scale` rather than the (now irrelevant) global
            # peak amplitude
            gained, gain = agc(raw, dt=tr.stats.delta, window_s=agc_window_s,
                               method=agc_method, eps_pct=agc_eps_pct,
                               return_gain=True)
            data = gained * wf_scale
            gain_db = 20 * np.log10(gain)
            gain_db_all.append(gain_db)
        else:
            data = (raw / global_max) * wf_scale
            gain_db = None

        times = tr.times()
        if shift_fn is not None:
            times = times - shift_fn(dist_km)

        processed[dist_km] = (times, data, gain_db)

    # Shared color scale for gain (dB), clipped to `gain_clip_pct` so a
    # handful of near-silent, extreme-gain samples don't wash out contrast
    # over the rest of the section
    norm, cmap = None, None
    if agc_window_s and gain_cmap:
        vmin, vmax = np.percentile(np.concatenate(gain_db_all), gain_clip_pct)
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap(gain_cmap)
        print(f"\tcoloring by AGC gain: [{vmin:.1f}, {vmax:.1f}]dB "
              f"(cmap='{gain_cmap}')")

    # Second pass: plot, using the shared gain color scale (if any)
    mappable = None
    for i, dist_km in enumerate(dists, 1):
        times, data, gain_db = processed[dist_km]
        y = data + dist_km

        if fill:
            fill_alpha = 0.2 if norm is not None else 0.4
            ax.fill_between(times, dist_km, y, where=(data >= 0),
                            color=wf_color, alpha=fill_alpha, zorder=4,
                            interpolate=True)

        if norm is not None:
            points = np.column_stack([times, y])[:, None, :]
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments, cmap=cmap, norm=norm, zorder=5)
            lc.set_array(gain_db[:-1])
            lc.set_linewidth(wf_lw)
            ax.add_collection(lc)
            mappable = lc
        else:
            ax.plot(times, y, c=wf_color, lw=wf_lw, zorder=5)

        print(f"\t[{i}/{len(dists)}] plotted d={dist_km:.1f}km")

    if mappable is not None:
        cbar = f.colorbar(mappable, ax=ax, pad=0.01)
        cbar.set_label("AGC gain applied [dB]")

    # Overlay TauP arrival-time curves so phase moveout/conversion is visible
    if tp_phases:
        cvals = cmaphex(nvals=len(tp_phases), cmap=tp_cmap)
        curves = taup_phase_curves(dists, src_depth_km, tp_phases,
                                   taup_model=taup_model,
                                   dist_step_km=dist_step_km, n_dist=n_dist)
        for i, phase in enumerate(tp_phases):
            dist_arr, time_arr = curves[phase]
            if len(dist_arr) == 0:
                continue
            if shift_fn is not None:
                # Shift the overlay curve consistently with the waveforms so
                # phases stay aligned to the correct reduced/aligned time
                time_arr = time_arr - shift_fn(dist_arr)
            ax.plot(time_arr, dist_arr, c=cvals[i], lw=tp_lw, ls="--",
                    zorder=10, label=phase)
            ax.text(time_arr[-1], dist_arr[-1], f" {phase}", c=cvals[i],
                    fontsize=9, va="center", zorder=11, fontweight="bold")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim(min(dists) - 5 * wf_scale, max(dists) + 5 * wf_scale)
    if title:
        ax.set_title(title)

    set_plot_aesthetic(ax, ytick_format="plain")

    # `set_plot_aesthetic` blanks y-ticks (prettyplot.py's YAXISOFF hack,
    # meant for momtenmeas.py's index-labeled plots), so re-apply explicit
    # distance labels here as a reference for which trace is which
    if ytick_step_km:
        dmin, dmax = min(dists), max(dists)
        yticks = np.arange(np.ceil(dmin / ytick_step_km) * ytick_step_km,
                           dmax + ytick_step_km, ytick_step_km)
    else:
        yticks = dists
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{d:g}" for d in yticks], fontsize=ytick_fontsize)

    plt.tight_layout()
    if save:
        print(f"saving to {save}")
        plt.savefig(save)
    if show:
        plt.show()
    plt.close(f)

    print(f"finished plotting in {time.time() - t0:.1f}s")

    return f, ax


def main(dist_km_list, baz=0, src_depth_km=10, tmin=1, tmax=4, corners=4,
         component="Z", syngine="iasp91_2s", taup_model="iasp91",
         path_cmtsolution=None, mt=None, tp_phases=None, wf_scale=8., xlim=None, ylim=None,
         reduce_by=None, agc_window_s=None, agc_method="rms",
         agc_eps_pct=1e-2, gain_cmap="inferno", gain_clip_pct=(1, 99),
         taper_vel_kmps=None, taper_vel_width_s=5.,
         ytick_step_km=None, fig_path="FIGURES", wav_path="SAC",
         parallel=True, show=True):
    """
    Gather Instaseis synthetics for a single source mechanism over a range of
    distances and plot the resulting distance record section
    """
    tp_phases = tp_phases or DEFAULT_PHASES
    t0 = time.time()
    print(f"=== DistRecSec: {len(dist_km_list)} distances; "
          f"syngine='{syngine}'; taup_model='{taup_model}'; "
          f"z={src_depth_km}km; baz={baz}deg; T=[{tmin}, {tmax}]s ===")

    drs = DistRecSec(dist_km_list=dist_km_list, baz=baz,
                     src_depth_km=src_depth_km, tmin=tmin, tmax=tmax,
                     corners=corners, components=component, syngine=syngine,
                     taup_model=taup_model, taper_vel_kmps=taper_vel_kmps,
                     taper_vel_width_s=taper_vel_width_s, fig_path=fig_path,
                     wav_path=wav_path)
    drs.get_src(path_cmtsolution=path_cmtsolution, mt=mt)

    save_tag = f"{syngine}_z{int(src_depth_km)}_b{int(baz % 360)}"
    dist_dict = drs.run(save_tag=save_tag, parallel=parallel)
    # `drs.run()` appends `drs.src_tag` (the mechanism identity) on top of
    # the context tag above, so read the actual resolved tag back here to
    # keep the figure filename in sync with the cache filenames
    save_tag = drs.save_tag

    title = (f"{drs.src_tag}  "
             f"{syngine}; z={src_depth_km}km; baz={int(baz % 360)}deg; "
             f"T=[{tmin}, {tmax}]s; comp={component}")
    if taper_vel_kmps:
        title += f"; tapered>{taper_vel_kmps}km/s"
    if reduce_by is not None:
        title += f"; reduced on {reduce_by}"
    save = f"{fig_path}/recsec_{save_tag}_{component}.png"

    plot_record_section(dist_dict, src_depth_km=src_depth_km,
                        tp_phases=tp_phases, taup_model=taup_model,
                        component=component, wf_scale=wf_scale, xlim=xlim,
                        ylim=ylim, reduce_by=reduce_by,
                        agc_window_s=agc_window_s,
                        agc_method=agc_method, agc_eps_pct=agc_eps_pct,
                        gain_cmap=gain_cmap, gain_clip_pct=gain_clip_pct,
                        ytick_step_km=ytick_step_km, title=title,
                        save=save, show=show)

    print(f"=== done in {time.time() - t0:.1f}s ===")


if __name__ == "__main__":
    main(
        dist_km_list=list(range(100, 2100, 100)),
        baz=0,
        src_depth_km=0.25,
        tmin=0.5, tmax=1, corners=8,
        xlim=[-10, 434],
        component="Z",
        syngine="ak135f_1s",
        taup_model="ak135f_no_mud",
        taper_vel_kmps=2.9,
        reduce_by="Pn",
        path_cmtsolution=f"CMTSOLUTION/CMTSOLUTION_{sys.argv[1]}",
        tp_phases=["Pn", "Pg", "Sn", "Sg"],
        wf_scale=30,
        agc_window_s=0,
        gain_cmap="viridis",
        parallel=True,
        show=True,
    )
