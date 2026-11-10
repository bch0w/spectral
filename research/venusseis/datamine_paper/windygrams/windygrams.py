#!/usr/bin/env python
"""
Plot ISHTAR (XX) waveforms, spectrograms and wind speed for high-, low- and
changing-wind windows at every site with an anemometer.

Each figure shows one component of one station (response removed,
velocity, highpassed), its spectrogram, and the 1-min mean wind speed from
the co-located anemometer (A) on a second axis. A window can be plotted for
all co-located stations of the site (XX Covered, Uncovered and Windscreen,
plus the permanent AK/HV station) and for the Z, N and E components.
Waveform y-limits, spectrogram colors and the wind axis are the same in
every figure so they can be compared directly; large transients are
clipped.

Everything that can be set (paths, sites, windows, what to plot,
processing and PrettyPlot options) is read from a YAML config file, by
default windygrams.yaml. A config can extend another with
``base: other.yaml`` and only list what it changes.

Figures are organized by the wind condition and window, e.g.
``FIGURES/high_wind/desd_high/XX.DESDU..HHZ.png``.

Windows are 4 h inside a single UTC day file, picked from 10-min wind
averages with find_windows.py (K24K and MITD hand-picked from the same
ranking). Wind noise dominates above ~5 Hz (15-60x louder in high wind);
below that, microseism noise hides it at MITD, hence the default highpass.

Example::

    python windygrams.py                      # windygrams.yaml
    python windygrams.py k24k_uz.yaml         # a smaller selection
    python windygrams.py k24k_uz.yaml --dry_run
"""
import argparse
import copy
import fnmatch
import os
import shlex
import sys
from zoneinfo import ZoneInfo

import matplotlib
import matplotlib.pyplot
import numpy as np
import yaml
from obspy import read, read_inventory, UTCDateTime

from pysep.pplot import PrettyPlot
from smooth_wind import smooth

FID_CONFIG = "windygrams.yaml"
STATIONS = ["C", "U", "W", "P"]
COMPONENTS = ["Z", "N", "E"]


def merge(base, new):
    """
    Recursively merge one config dictionary onto another. Nested
    dictionaries are merged key by key; anything else in `new` (including
    lists) replaces the value in `base`

    :type base: dict
    :param base: configuration to start from, not modified
    :type new: dict
    :param new: configuration whose values take precedence
    :rtype: dict
    :return: merged configuration
    """
    out = copy.deepcopy(base)
    for key, val in new.items():
        if isinstance(val, dict) and isinstance(out.get(key), dict):
            out[key] = merge(out[key], val)
        else:
            out[key] = copy.deepcopy(val)
    return out


def read_config(fid):
    """
    Read a YAML config file, following its `base` chain if it has one

    :type fid: str
    :param fid: path to the YAML file
    :rtype: dict
    :return: configuration with all bases merged in
    """
    with open(fid) as f:
        cfg = yaml.safe_load(f) or {}
    base = cfg.pop("base", None)
    if base:
        # A relative base path is relative to the file that names it
        fid_base = os.path.join(os.path.dirname(os.path.abspath(fid)), base)
        cfg = merge(read_config(fid_base), cfg)
    return cfg


def check_config(cfg):
    """
    Check a configuration and normalize its values in place: numbers are
    converted to float (YAML reads e.g. '2E-5' as a string), window times
    must be quoted strings, and window and station choices must exist

    :type cfg: dict
    :param cfg: configuration from `read_config`
    :raises ValueError: for a missing section or an invalid value
    """
    for key in ["paths", "year", "select", "processing", "pplot", "sites",
                "windows"]:
        if key not in cfg:
            raise ValueError(f"config is missing the '{key}' section")

    proc = cfg["processing"]
    for key in ["fmin", "ylim", "wind_smooth", "env_smooth"]:
        proc[key] = float(proc[key])
    proc["envelope"] = bool(proc.get("envelope", False))

    for name, win in cfg["windows"].items():
        for key in ["site", "jday", "start", "end", "desc"]:
            if key not in win:
                raise ValueError(f"window '{name}' is missing '{key}'")
        if win["site"] not in cfg["sites"]:
            raise ValueError(f"window '{name}': unknown site "
                             f"'{win['site']}'")
        for key in ["start", "end"]:
            if not isinstance(win[key], str):
                raise ValueError(f"window '{name}': quote the {key} time, "
                                 f"e.g. \"13:10\"; YAML read it as "
                                 f"{win[key]!r}")

    sel = cfg["select"]
    for key, choices in [("stations", STATIONS),
                         ("components", COMPONENTS)]:
        bad = set(sel[key]) - set(choices)
        if bad:
            raise ValueError(f"select.{key}: unknown {sorted(bad)}, "
                             f"choose from {choices}")


def select_windows(cfg):
    """
    Return the window names chosen by `select.windows`, which may be 'all',
    a name, or a list of names and wildcards such as 'k24k_*'

    :type cfg: dict
    :param cfg: configuration
    :rtype: list of str
    :return: window names in the order of the config
    :raises ValueError: if a name or wildcard matches no window
    """
    patterns = cfg["select"]["windows"]
    if patterns == "all":
        return list(cfg["windows"])
    if isinstance(patterns, str):
        patterns = [patterns]
    names = []
    for pattern in patterns:
        match = fnmatch.filter(cfg["windows"], pattern)
        if not match:
            raise ValueError(f"select.windows: '{pattern}' matches no "
                             f"window")
        names += [m for m in match if m not in names]
    return names


def wind_file(cfg, site, jday):
    """
    Return the smoothed wind speed file for a site and day, making it from
    the raw 1 Hz anemometer file if it does not exist yet

    :type cfg: dict
    :param cfg: configuration
    :type site: str
    :param site: site name, e.g. 'K24K'
    :type jday: int
    :param jday: julian day
    :rtype: str
    :return: path to the smoothed wind speed file
    """
    name = f"XX.{site}A..LWS.{cfg['year']}.{jday:03d}"
    fid = os.path.join(cfg["paths"]["wind"], name)
    if not os.path.exists(fid):
        smooth(os.path.join(cfg["paths"]["data"], name),
               cfg["processing"]["wind_smooth"], cfg["paths"]["wind"])
    return fid


def title(cfg, code, win, fid_wind, start, end):
    """
    Build a figure title with the station, wind summary, local time window
    and filter band

    :type cfg: dict
    :param cfg: configuration
    :type code: str
    :param code: seismic trace ID, e.g. 'XX.K24KU..EHZ'
    :type win: dict
    :param win: window entry of the config
    :type fid_wind: str
    :param fid_wind: smoothed wind speed file
    :type start: obspy.UTCDateTime
    :param start: window start
    :type end: obspy.UTCDateTime
    :param end: window end
    :rtype: str
    :return: title text
    """
    wind = read(fid_wind, starttime=start, endtime=end)[0].data
    tz = ZoneInfo(cfg["sites"][win["site"]]["tz"])
    t0 = start.datetime.replace(tzinfo=ZoneInfo("UTC")).astimezone(tz)
    t1 = end.datetime.replace(tzinfo=ZoneInfo("UTC")).astimezone(tz)
    fmin = cfg["processing"]["fmin"]
    # The temporary XX network code is left out of figure labels
    label = code.removeprefix("XX.")
    return (f"{label} {win['desc']} (wind {np.mean(wind):.1f} m/s mean), "
            f"{t0:%Y-%m-%d %H:%M}-{t1:%H:%M} {t0:%Z}, >{fmin:g} Hz")


def trace_id(cfg, site, sta, comp):
    """
    Return the trace ID of one component of a co-located station

    :type cfg: dict
    :param cfg: configuration
    :type site: str
    :param site: site name, key of `sites` in the config
    :type sta: str
    :param sta: C, U or W for the XX stations, P for the permanent station
    :type comp: str
    :param comp: component, Z, N or E
    :rtype: str
    :return: trace ID, e.g. 'XX.DESDU..HHZ' or 'HV.DESD..EHZ'
    """
    info = cfg["sites"][site]
    if sta == "P":
        return f"{info['perm']}..{info['perm_cha']}{comp}"
    return f"XX.{site}{sta}..{info['cha']}{comp}"


def kind_dir(desc):
    """
    Return the output subdirectory for a window description

    :type desc: str
    :param desc: window description, e.g. 'high wind'
    :rtype: str
    :return: directory name, e.g. 'high_wind'
    """
    return desc.replace(" ", "_")


def read_inventories(cfg):
    """
    Read the ISHTAR StationXML and that of every permanent station

    :type cfg: dict
    :param cfg: configuration
    :rtype: obspy.core.inventory.Inventory
    :return: one inventory with all stations
    """
    path = cfg["paths"]["xml"]
    inv = read_inventory(os.path.join(path, "ISHTAR.xml"))
    for perm in sorted({site["perm"] for site in cfg["sites"].values()}):
        inv += read_inventory(os.path.join(path, f"{perm}.xml"))
    return inv


def figure_setup(cfg, name, sta="U", comp="Z"):
    """
    Work out the input files, output file and PrettyPlot keyword arguments
    of one figure, without reading any waveform data

    :type cfg: dict
    :param cfg: configuration
    :type name: str
    :param name: window name, key of `windows` in the config
    :type sta: str
    :param sta: station: C (covered), U (uncovered), W (windscreen) or P
        (permanent AK/HV station)
    :type comp: str
    :param comp: component, Z, N or E
    :rtype: dict or None
    :return: 'code', 'fid_data', 'fid_wind', 'fid_raw_wind', 'fid_inv',
        'other_cha' (other channels in the day file) and 'kwargs' for
        PrettyPlot (without the inventory); None if there is no data file
    """
    win = cfg["windows"][name]
    proc = cfg["processing"]
    year, jday = cfg["year"], int(win["jday"])
    site = win["site"]
    code = trace_id(cfg, site, sta, comp)
    cha = code.split(".")[-1]
    fid_data = os.path.join(cfg["paths"]["data"],
                            f"{code}.{year}.{jday:03d}")
    if not os.path.exists(fid_data):
        return None
    # Header only: sampling rate for the pre-filter, other channels
    st = read(fid_data, headonly=True)
    samp_rate = st.select(channel=cha)[0].stats.sampling_rate
    other_cha = sorted({tr.stats.channel for tr in st} - {cha})

    start = UTCDateTime(f"{year}-{jday:03d}T{win['start']}")
    end = UTCDateTime(f"{year}-{jday:03d}T{win['end']}")
    wind = f"XX.{site}A..LWS.{year}.{jday:03d}"
    fid_wind = os.path.join(cfg["paths"]["wind"], wind)
    fid_raw_wind = os.path.join(cfg["paths"]["data"], wind)
    perm = cfg["sites"][site]["perm"]
    fid_inv = os.path.join(cfg["paths"]["xml"],
                           f"{perm}.xml" if sta == "P" else "ISHTAR.xml")
    path_fig = os.path.join(cfg["paths"]["out"], kind_dir(win["desc"]),
                            name)

    kwargs = dict(cfg["pplot"])
    suffix = ""
    # The envelope is drawn as a solid line over the normal waveform
    if proc["envelope"]:
        kwargs.update(wf_env=True, wf_env_smooth=proc["env_smooth"])
        suffix = "_env"
    else:
        kwargs.pop("wf_env_color", None)
    nyq = samp_rate / 2
    kwargs.update(
        fmin=proc["fmin"], ylim=[-proc["ylim"], proc["ylim"]],
        xlim=[str(start), str(end)],
        # Pre-filter the response removal up to Nyquist of this channel
        rr_pre_filt=[0.05, 0.1, 0.9 * nyq, nyq],
        labels=[code.removeprefix("XX.")],
        add_trace=fid_wind, tr_label=f"{site}A wind speed",
        title=title(cfg, code, win, fid_wind if os.path.exists(fid_wind)
                    else fid_raw_wind, start, end),
        save=os.path.join(path_fig, f"{code}{suffix}.png"),
    )
    return {"code": code, "fid_data": fid_data, "fid_wind": fid_wind,
            "fid_raw_wind": fid_raw_wind, "fid_inv": fid_inv,
            "other_cha": other_cha, "kwargs": kwargs}


# PrettyPlot keyword arguments whose command line flag is the negation
NEGATED_FLAGS = {"show": "--noshow", "legend": "--no_legend",
                 "sp_dbscale": "--sp_dbscale_off",
                 "rr_zero_mean": "--rr_zero_mean_off",
                 "rr_taper": "--rr_taper_off"}


def _cli_value(val):
    """
    Format one value for the pplot command line. Floats are written in
    plain decimal so argparse does not take e.g. '-2e-05' for an option

    :type val: any
    :param val: value of a PrettyPlot keyword argument
    :rtype: str
    :return: shell-quoted text
    """
    if isinstance(val, float):
        val = np.format_float_positional(val, trim="-")
    return shlex.quote(str(val))


def pplot_command(setup, wind_smooth=60., python=None):
    """
    Build the shell command(s) that make one figure with the pplot command
    line tool, equivalent to `plot_window`

    :type setup: dict
    :param setup: output of `figure_setup`
    :type wind_smooth: float
    :param wind_smooth: running mean of the wind speed in seconds, used if
        the smoothed wind file still has to be made
    :type python: str or None
    :param python: interpreter to run smooth_wind.py with, default the one
        running this script
    :rtype: str
    :return: one command, preceded by the smooth_wind.py call if the
        smoothed wind file does not exist yet
    """
    lines = []
    if not os.path.exists(setup["fid_wind"]):
        python = python or sys.executable
        lines.append(f"{python} smooth_wind.py {setup['fid_raw_wind']} "
                     f"--window {wind_smooth:g} "
                     f"--path_out {os.path.dirname(setup['fid_wind'])}")
    if setup["other_cha"]:
        lines.append(f"# {setup['fid_data']} also holds "
                     f"{', '.join(setup['other_cha'])}, which pplot will "
                     f"plot too; windygrams.py drops them")
    args = [setup["fid_data"], "--inv", setup["fid_inv"]]
    for key, val in setup["kwargs"].items():
        if key in NEGATED_FLAGS:
            if not val:
                args.append(NEGATED_FLAGS[key])
        elif isinstance(val, bool):
            if val:
                args.append(f"--{key}")
        elif val is None:
            continue
        elif isinstance(val, (list, tuple)):
            args += [f"--{key}"] + [_cli_value(v) for v in val]
        else:
            args += [f"--{key}", _cli_value(val)]
    lines.append("pplot " + " ".join(args))
    return "\n".join(lines)


def plot_window(cfg, name, inv, sta="U", comp="Z"):
    """
    Plot one window for one component of one station with PrettyPlot

    :type cfg: dict
    :param cfg: configuration
    :type name: str
    :param name: window name, key of `windows` in the config
    :type inv: obspy.core.inventory.Inventory
    :param inv: StationXML inventory used for response removal
    :type sta: str
    :param sta: station: C (covered), U (uncovered), W (windscreen) or P
        (permanent AK/HV station)
    :type comp: str
    :param comp: component, Z, N or E
    :rtype: str or None
    :return: path of the saved figure, or None if there is no data
    """
    code = trace_id(cfg, cfg["windows"][name]["site"], sta, comp)
    setup = figure_setup(cfg, name, sta, comp)
    if setup is None:
        print(f"   no data for {code}")
        return None
    fid_data, kwargs = setup["fid_data"], setup["kwargs"]
    # Make the smoothed wind file before plotting if needed
    wind_file(cfg, cfg["windows"][name]["site"],
              int(cfg["windows"][name]["jday"]))
    os.makedirs(os.path.dirname(kwargs["save"]), exist_ok=True)

    pp = PrettyPlot([fid_data], inv=inv, **kwargs)
    # Some day files hold other channels too (e.g. AK LHZ next to BHZ)
    pp.st = pp.st.select(channel=code.split(".")[-1])
    try:
        pp.run()
    except SystemExit:
        # PrettyPlot exits if trimming leaves no data in the window
        print(f"   no data in window: {fid_data}")
        return None
    except Exception as e:
        print(f"   failed: {fid_data}: {e}")
        matplotlib.pyplot.close("all")
        return None
    return kwargs["save"]


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", nargs="?", default=FID_CONFIG,
                        help=f"YAML config file, default {FID_CONFIG}")
    parser.add_argument("--dry_run", action="store_true",
                        help="list the figures the config would make, "
                             "with the pplot command for each, without "
                             "plotting")
    return parser.parse_args()


def main():
    """Plot every selected window, station and component of a config"""
    args = parse_args()
    fid_config = os.path.abspath(args.config)
    # Paths in the config are relative to this script's directory
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    cfg = read_config(fid_config)
    check_config(cfg)
    names = select_windows(cfg)
    sel = cfg["select"]
    print(f"{args.config}: {len(names)} window(s) x "
          f"{len(sel['stations'])} station(s) x "
          f"{len(sel['components'])} component(s)")

    if args.dry_run:
        print(f"# pplot commands, to run from {os.getcwd()}")
        for name in names:
            win = cfg["windows"][name]
            print(f"\n# {name} ({win['desc']})")
            for sta in sel["stations"]:
                for comp in sel["components"]:
                    setup = figure_setup(cfg, name, sta, comp)
                    if setup is None:
                        code = trace_id(cfg, win["site"], sta, comp)
                        print(f"# no data for {code}")
                    else:
                        print(pplot_command(
                            setup, cfg["processing"]["wind_smooth"]))
        return

    matplotlib.use("Agg")
    # Hours of 100 Hz data exceed Agg's single-path limit otherwise
    matplotlib.rcParams["agg.path.chunksize"] = 10000
    inv = read_inventories(cfg)
    for name in names:
        for sta in sel["stations"]:
            for comp in sel["components"]:
                fid = plot_window(cfg, name, inv, sta=sta, comp=comp)
                if fid:
                    print(f"-> {fid}")


if __name__ == "__main__":
    main()
