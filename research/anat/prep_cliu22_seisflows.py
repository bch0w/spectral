#!/usr/bin/env python3
"""
Prepare the Liu et al. (2022) empirical Green's function dataset for SeisFlows

Converts the Zenodo collection of Liu et al. (2022) Empirical Green's
Functions (three-station ambient noise cross-correlations) from its native
format into a directory structure and SAC header convention that SeisFlows
can use. The original data (the directory extracted from the Zenodo .zip)
are only ever READ. Every modified file is written into a new output directory,
so the original dataset is never touched.

The workflow is split into subcommands:

1. ``inspect``: print the SAC headers of a few input files so the time
   reference (zero-lag time, ``b``, ``o``, one- or two-sided) can be checked
   BEFORE anything is written
2. ``prep``: copy one data type ('hyp' by default) into the output
   directory, optionally keeping only the sources and receivers listed in
   SPECFEM STATIONS files. Each file is renamed, its SAC headers are fixed,
   it is shifted to a new origin time and trimmed to a fixed length. The
   output is structured as

       <OUTPUT_DIR>/<NET>_<STA>/<ZZ|TT>/<NET>.<STA>.LX<Z|T>.SAC

   where the directory names the (virtual) source station and the file
   names the receiver station
3. ``select``: count measurements for each source and receiver in a prepped
   directory, write thresholded SOURCES and STATIONS files, and report
   station codes that are shared by more than one network (e.g., AK and TA)
4. ``plot``: make station coverage figures for a prepped directory

Noteworthy points:

- 'hyp' (hyperbolic) data are used by default because they have higher SNR
  than 'ell' (elliptical) data (Liu et al. 2022)
- Channels are named 'LX?': L because the data are sampled at 1 Hz, X
  because the time series is derived from observational data
- Origin time is kept at the zero-lag time of the Zenodo data
  (2019-01-01T00:00:00) by default. It can be shifted by setting
  ``--new-origin`` to something other than ``--old-origin``. Check the
  zero-lag time with ``inspect`` before running ``prep``

Example::

    python prep_cliu22_seisflows.py inspect SAC_I3_stack_4_Zendo
    python prep_cliu22_seisflows.py prep SAC_I3_stack_4_Zendo LIU22_EGF \\
        --sources SOURCES_NALASKA --stations STATIONS_NALASKA
    python prep_cliu22_seisflows.py select LIU22_EGF STATIONS_ALL \\
        --src-threshold 86 --rcv-threshold 5
    python prep_cliu22_seisflows.py plot LIU22_EGF STATIONS_ALL
"""
import argparse
import os
import sys
from collections import Counter, defaultdict
from glob import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from obspy import read, UTCDateTime  # noqa: E402
from obspy.core import AttribDict  # noqa: E402
from obspy.io.sac.util import get_sac_reftime  # noqa: E402
from pysep.utils.io import read_stations, write_stations_file  # noqa: E402


# Stack directories are named e.g., 'Lov_I3_hyp_stack'
KERNELS = {"Lov": "TT", "Ray": "ZZ"}
OLD_ORIGIN = "2019-01-01T00:00:00"
NEW_ORIGIN = "2019-01-01T00:00:00"
LENGTH_S = 20 * 60.


# =============================================================================
#                              READING ORIGINAL DATA
# =============================================================================
def iter_input_files(input_dir):
    """
    Yield every data file in the original (extracted) dataset

    Data files are expected at ``<INPUT_DIR>/<STACK>/<SRC>/<FILE>``, e.g.,
    ``SAC_I3_stack_4_Zendo/Lov_I3_hyp_stack/AK_ANM/<FILE>.SAC``. Files are
    only listed here, never modified.

    :type input_dir: str
    :param input_dir: extracted Zenodo directory, e.g.,
        'SAC_I3_stack_4_Zendo'
    :rtype: generator of str
    :return: path to each data file, sorted alphabetically
    :raises FileNotFoundError: if ``input_dir`` is not a directory
    """
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"{input_dir} is not a directory")
    for path in sorted(glob(os.path.join(input_dir, "*", "*", "*"))):
        if os.path.isfile(path):
            yield path


def split_input_path(path):
    """
    Split a data file path into its stack, source and filename components

    :type path: str
    :param path: path to a data file from :func:`iter_input_files`
    :rtype: tuple of str
    :return: (stack, source, filename), e.g., ('Lov_I3_hyp_stack',
        'AK_ANM', '<FILE>.SAC')
    """
    src_dir, fname = os.path.split(path)
    stack_dir, src = os.path.split(src_dir)
    return os.path.basename(stack_dir), src, fname


def parse_stack_name(stack):
    """
    Get the kernel, component and data type from a stack directory name

    :type stack: str
    :param stack: stack directory name, e.g., 'Lov_I3_hyp_stack'
    :rtype: tuple of str or None
    :return: (kernel, component, data type), e.g., ('TT', 'T', 'hyp'), or
        None if the name is not a recognized stack directory
    """
    parts = stack.split("_")
    if len(parts) != 4 or parts[0] not in KERNELS:
        return None
    kernel = KERNELS[parts[0]]
    return kernel, kernel[-1], parts[2].lower()


def parse_data_filename(fname):
    """
    Get the source and receiver station codes from a data filename

    Filenames are formatted ``<PREFIX>_<NETSRC>_<STASRC>_<NETRCV>_<STARCV>``
    followed by a file extension.

    :type fname: str
    :param fname: basename of the data file
    :rtype: tuple of str
    :return: (source network, source station, receiver network, receiver
        station)
    :raises ValueError: if the filename does not have the expected format
    """
    stem = os.path.splitext(fname)[0]
    parts = stem.split("_")
    if len(parts) != 5:
        raise ValueError(f"unexpected filename format: {fname}")
    return tuple(parts[1:])


def read_station_codes(stations_file):
    """
    Read a SPECFEM STATIONS file into a set of (network, station) codes

    :type stations_file: str
    :param stations_file: path to a SPECFEM-formatted STATIONS file
    :rtype: set of tuple of str
    :return: (network, station) codes listed in the file
    """
    inv = read_stations(stations_file)
    return {(net.code, sta.code) for net in inv for sta in net}


def read_station_coords(stations_file):
    """
    Read a SPECFEM STATIONS file into a coordinate lookup table

    :type stations_file: str
    :param stations_file: path to a SPECFEM-formatted STATIONS file
    :rtype: dict
    :return: {'NET_STA': (latitude, longitude)} for each station in the file
    """
    inv = read_stations(stations_file)
    return {f"{net.code}_{sta.code}": (sta.latitude, sta.longitude)
            for net in inv for sta in net}


# =============================================================================
#                                    INSPECT
# =============================================================================
def inspect(input_dir, nfiles=3, old_origin=OLD_ORIGIN):
    """
    Print SAC headers of the first few files of each stack directory

    Used to confirm the time reference of the original data before running
    :func:`prep`. Nothing is written.

    :type input_dir: str
    :param input_dir: extracted Zenodo directory
    :type nfiles: int
    :param nfiles: number of files to print for each stack directory
    :type old_origin: str
    :param old_origin: expected zero-lag time of the original data
    """
    old_origin = UTCDateTime(old_origin)
    shown = Counter()
    for path in iter_input_files(input_dir):
        stack, _, _ = split_input_path(path)
        if shown[stack] >= nfiles:
            continue
        shown[stack] += 1

        tr = read(path, format="SAC")[0]
        sac = tr.stats.sac
        print(f"\n{path}")
        print(f"    npts={tr.stats.npts}  delta={tr.stats.delta}")
        print(f"    starttime={tr.stats.starttime}  "
              f"endtime={tr.stats.endtime}")
        print(f"    reftime={get_sac_reftime(sac)}  b={sac.get('b')}  "
              f"e={sac.get('e')}  o={sac.get('o')}")
        print(f"    knetwk={sac.get('knetwk')}  kstnm={sac.get('kstnm')}  "
              f"kcmpnm={sac.get('kcmpnm')}")
        print(f"    evla={sac.get('evla')}  evlo={sac.get('evlo')}  "
              f"stla={sac.get('stla')}  stlo={sac.get('stlo')}  "
              f"dist={sac.get('dist')}")
        lag0 = tr.stats.starttime - old_origin
        side = "TWO-SIDED" if lag0 < 0 else "one-sided"
        print(f"    first sample is at lag {lag0:.2f}s w.r.t. old origin "
              f"{old_origin} -> {side}")

    if not shown:
        print(f"No data files found in {input_dir}")


# =============================================================================
#                                      PREP
# =============================================================================
def fix_trace(tr, network, station, channel, old_origin, new_origin,
              length_s):
    """
    Rename a trace, shift it to a new origin time and trim it

    The trace is shifted so that the sample at ``old_origin`` (zero lag) ends
    up at ``new_origin``, then trimmed to [new_origin, new_origin +
    length_s], i.e., only positive lags are kept. The SAC reference time is
    set to ``new_origin`` so that ``b`` and ``o`` are both relative to it.

    :type tr: obspy.core.trace.Trace
    :param tr: trace to modify in place
    :type network: str
    :param network: new network code
    :type station: str
    :param station: new station code
    :type channel: str
    :param channel: new channel code
    :type old_origin: obspy.UTCDateTime
    :param old_origin: zero-lag time of the original data
    :type new_origin: obspy.UTCDateTime
    :param new_origin: zero-lag time of the output data
    :type length_s: float
    :param length_s: length of the output trace in seconds
    :raises ValueError: if the trace starts after ``old_origin``, i.e., the
        zero-lag sample is missing and ``old_origin`` is probably wrong
    """
    if tr.stats.starttime - old_origin > tr.stats.delta / 2:
        raise ValueError(f"trace starts at {tr.stats.starttime}, after the "
                         f"old origin {old_origin}")

    tr.stats.network = network
    tr.stats.station = station
    tr.stats.channel = channel

    tr.stats.starttime += new_origin - old_origin
    tr.trim(new_origin, new_origin + length_s, nearest_sample=True)

    if "sac" not in tr.stats:
        tr.stats.sac = AttribDict()
    sac = tr.stats.sac
    sac.nzyear = new_origin.year
    sac.nzjday = new_origin.julday
    sac.nzhour = new_origin.hour
    sac.nzmin = new_origin.minute
    sac.nzsec = new_origin.second
    sac.nzmsec = new_origin.microsecond // 1000
    sac.o = 0.
    sac.evdp = 0.


def _write_atomic(st, fid_out):
    """
    Write a Stream as SAC via a temporary file so partial files never exist

    :type st: obspy.core.stream.Stream
    :param st: stream to write
    :type fid_out: str
    :param fid_out: final output path
    """
    os.makedirs(os.path.dirname(fid_out), exist_ok=True)
    fid_tmp = f"{fid_out}.tmp"
    st.write(fid_tmp, format="SAC")
    os.replace(fid_tmp, fid_out)


def _check_output_dir(input_dir, output_dir):
    """
    Make sure the output directory cannot overwrite the original data

    :type input_dir: str
    :param input_dir: extracted Zenodo directory
    :type output_dir: str
    :param output_dir: directory that modified files will be written to
    :raises ValueError: if the output directory is, or is inside, the
        original data directory
    """
    input_abs = os.path.realpath(input_dir)
    output_abs = os.path.realpath(output_dir)
    if output_abs == input_abs or output_abs.startswith(input_abs + os.sep):
        raise ValueError("output directory must not be inside the input "
                         "directory")


def prep(input_dir, output_dir, data_type="hyp", sources_file=None,
         stations_file=None, old_origin=OLD_ORIGIN, new_origin=NEW_ORIGIN,
         length_s=LENGTH_S, overwrite=False, dry_run=False):
    """
    Copy, rename, fix headers of, shift and trim the original dataset

    Original files are only read and the modified copies are written to
    ``output_dir``. Writes are atomic, so a file that exists in the output
    directory has been fully processed and re-running only processes files
    that are missing (or all files if ``overwrite``). Failures are listed in
    'prep_failures.txt' in the output directory.

    :type input_dir: str
    :param input_dir: extracted Zenodo directory
    :type output_dir: str
    :param output_dir: directory to write the modified dataset to
    :type data_type: str
    :param data_type: 'hyp' or 'ell'
    :type sources_file: str or None
    :param sources_file: optional STATIONS file; source stations not listed
        are skipped
    :type stations_file: str or None
    :param stations_file: optional STATIONS file; receiver stations not
        listed are skipped
    :type old_origin: str
    :param old_origin: zero-lag time of the original data
    :type new_origin: str
    :param new_origin: zero-lag time of the output data
    :type length_s: float
    :param length_s: length of output traces in seconds
    :type overwrite: bool
    :param overwrite: re-process files that already exist in the output
    :type dry_run: bool
    :param dry_run: only count what would be written, write nothing
    """
    _check_output_dir(input_dir, output_dir)
    old_origin = UTCDateTime(old_origin)
    new_origin = UTCDateTime(new_origin)
    src_codes = read_station_codes(sources_file) if sources_file else None
    rcv_codes = read_station_codes(stations_file) if stations_file else None

    counts = Counter()
    failures = []
    for path in iter_input_files(input_dir):
        stack, src_dir, fname = split_input_path(path)
        stack_info = parse_stack_name(stack)
        if stack_info is None:
            counts["unrecognized stack dir"] += 1
            continue
        kernel, comp, type_ = stack_info
        if type_ != data_type:
            counts[f"skipped type '{type_}'"] += 1
            continue

        try:
            net_src, sta_src, net_rcv, sta_rcv = parse_data_filename(fname)
        except ValueError as e:
            failures.append((path, str(e)))
            continue
        if src_dir != f"{net_src}_{sta_src}":
            counts["source dir/filename mismatch (used filename)"] += 1

        if src_codes is not None and (net_src, sta_src) not in src_codes:
            counts["skipped source not in sources file"] += 1
            continue
        if rcv_codes is not None and (net_rcv, sta_rcv) not in rcv_codes:
            counts["skipped receiver not in stations file"] += 1
            continue

        channel = f"LX{comp}"
        fid_out = os.path.join(output_dir, f"{net_src}_{sta_src}", kernel,
                               f"{net_rcv}.{sta_rcv}.{channel}.SAC")
        if os.path.exists(fid_out) and not overwrite:
            counts["already exists"] += 1
            continue
        if dry_run:
            counts["would write"] += 1
            continue

        try:
            st = read(path, format="SAC")
            if len(st) != 1:
                raise ValueError(f"expected 1 trace, found {len(st)}")
            fix_trace(st[0], net_rcv, sta_rcv, channel, old_origin,
                      new_origin, length_s)
            _write_atomic(st, fid_out)
        except Exception as e:
            failures.append((path, str(e)))
            continue
        counts["written"] += 1
        if counts["written"] % 1000 == 0:
            print(f"{counts['written']} files written")

    print(f"\nSummary for {input_dir} -> {output_dir}")
    for key, val in sorted(counts.items()):
        print(f"    {key}: {val}")
    print(f"    failed: {len(failures)}")
    if failures:
        os.makedirs(output_dir, exist_ok=True)
        fid = os.path.join(output_dir, "prep_failures.txt")
        with open(fid, "w") as f:
            for path, err in failures:
                f.write(f"{path}: {err}\n")
        print(f"    failures written to {fid}")


# =============================================================================
#                                STATION SELECTION
# =============================================================================
def iter_prepped_files(data_dir):
    """
    Yield each source, kernel and receiver in a prepped directory

    :type data_dir: str
    :param data_dir: output directory of :func:`prep`
    :rtype: generator of tuple of str
    :return: ('NET_STA' source code, kernel, 'NET_STA' receiver code)
    """
    for path in sorted(glob(os.path.join(data_dir, "*_*", "*", "*.SAC"))):
        kernel_dir, fname = os.path.split(path)
        src_dir, kernel = os.path.split(kernel_dir)
        rcv_net, rcv_sta, *_ = fname.split(".")
        yield os.path.basename(src_dir), kernel, f"{rcv_net}_{rcv_sta}"


def count_source_receiver_hits(data_dir, src_threshold=86, rcv_threshold=5):
    """
    Count the number of measurements for each source and receiver station

    Used to decide which stations to keep. Stations with fewer than the
    threshold number of measurements are excluded from the returned counts.

    :type data_dir: str
    :param data_dir: output directory of :func:`prep`
    :type src_threshold: int
    :param src_threshold: minimum number of measurements per source, 0 to
        keep all
    :type rcv_threshold: int
    :param rcv_threshold: minimum number of measurements per receiver, 0 to
        keep all
    :rtype: tuple of dict
    :return: ({'NET_STA': count} for sources, {'NET_STA': count} for
        receivers), each sorted by increasing count
    """
    src_count, rcv_count = Counter(), Counter()
    for src_code, _, rcv_code in iter_prepped_files(data_dir):
        src_count[src_code] += 1
        rcv_count[rcv_code] += 1

    src_count = {k: v for k, v in sorted(src_count.items(),
                                         key=lambda item: item[1])
                 if v >= src_threshold}
    rcv_count = {k: v for k, v in sorted(rcv_count.items(),
                                         key=lambda item: item[1])
                 if v >= rcv_threshold}
    return src_count, rcv_count


def find_doubled_sources(src_count):
    """
    Find source station codes that are shared by more than one network

    AK and TA stations share locations but have different datasets. One of
    each pair should be removed to avoid running two simulations with
    almost exactly the same source.

    :type src_count: dict
    :param src_count: {'NET_STA': count} from
        :func:`count_source_receiver_hits`
    :rtype: dict
    :return: {'STA': [('NET_STA', count), ...]} sorted by station code, with
        each list sorted by decreasing count
    """
    by_sta = defaultdict(list)
    for code, count in src_count.items():
        by_sta[code.split("_")[1]].append((code, count))
    return {sta: sorted(codes, key=lambda item: -item[1])
            for sta, codes in sorted(by_sta.items()) if len(codes) > 1}


def select_stations(data_dir, stations_all, src_threshold=86,
                    rcv_threshold=5, output_dir="."):
    """
    Write thresholded SOURCES and STATIONS files and a doubled-source report

    :type data_dir: str
    :param data_dir: output directory of :func:`prep`
    :type stations_all: str
    :param stations_all: STATIONS file with coordinates of every station
    :type src_threshold: int
    :param src_threshold: minimum number of measurements per source
    :type rcv_threshold: int
    :param rcv_threshold: minimum number of measurements per receiver
    :type output_dir: str
    :param output_dir: directory to write 'SOURCES', 'STATIONS' and
        'check_doubles.txt' to
    """
    src_count, rcv_count = count_source_receiver_hits(
        data_dir, src_threshold, rcv_threshold)
    inv = read_stations(stations_all)

    src_inv, rcv_inv = inv.copy(), inv.copy()
    for net in inv:
        for sta in net:
            code = f"{net.code}_{sta.code}"
            if code not in src_count:
                src_inv = src_inv.remove(network=net.code, station=sta.code)
            if code not in rcv_count:
                rcv_inv = rcv_inv.remove(network=net.code, station=sta.code)

    missing = (set(src_count) | set(rcv_count)) - {
        f"{net.code}_{sta.code}" for net in inv for sta in net}
    if missing:
        print(f"WARNING: {len(missing)} stations with data are not in "
              f"{stations_all}: {sorted(missing)}")

    print(f"Sources: {len(src_count)} (>= {src_threshold} measurements)")
    print(f"Receivers: {len(rcv_count)} (>= {rcv_threshold} measurements)")

    os.makedirs(output_dir, exist_ok=True)
    write_stations_file(src_inv, os.path.join(output_dir, "SOURCES"))
    write_stations_file(rcv_inv, os.path.join(output_dir, "STATIONS"))

    doubles = find_doubled_sources(src_count)
    fid = os.path.join(output_dir, "check_doubles.txt")
    with open(fid, "w") as f:
        for sta, codes in doubles.items():
            line = ", ".join(f"{code}={count}" for code, count in codes)
            f.write(f"{sta}: {line}\n")
    print(f"{len(doubles)} doubled source station codes written to {fid}")


# =============================================================================
#                                   PLOTTING
# =============================================================================
def _set_map_axes(ax, extent=None):
    """
    Label axes of a longitude/latitude map and optionally set its extent

    :type ax: matplotlib.axes.Axes
    :param ax: axis to modify
    :type extent: list of float or None
    :param extent: [lon_min, lon_max, lat_min, lat_max], None for automatic
    """
    if extent is not None:
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")


def plot_source_receiver_hits(src_count, rcv_count, coords, output_dir,
                              extent=None):
    """
    Map the number of measurements at each source and receiver station

    Stations in ``coords`` that are not counted (e.g., thresholded out) are
    plotted as gray crosses for reference.

    :type src_count: dict
    :param src_count: {'NET_STA': count} for sources
    :type rcv_count: dict
    :param rcv_count: {'NET_STA': count} for receivers
    :type coords: dict
    :param coords: {'NET_STA': (lat, lon)} from :func:`read_station_coords`
    :type output_dir: str
    :param output_dir: directory to save figures to
    :type extent: list of float or None
    :param extent: [lon_min, lon_max, lat_min, lat_max], None for automatic
    """
    for key, counts in [("source", src_count), ("receiver", rcv_count)]:
        codes = [code for code in counts if code in coords]
        if not codes:
            print(f"No {key} stations to plot")
            continue
        lats, lons = zip(*[coords[code] for code in codes])
        cnts = [counts[code] for code in codes]

        f, ax = plt.subplots(figsize=(20, 10), dpi=200)
        sc = ax.scatter(lons, lats, c=cnts, marker="o", zorder=6, s=40,
                        ec="k", lw=1)
        for lon, lat, code in zip(lons, lats, codes):
            ax.text(s=f"{code}: {counts[code]}", x=lon, y=lat, size=7,
                    zorder=8)
        f.colorbar(sc, ax=ax, label="counts")

        # Plot all stations for reference of what has been kicked
        for code, (lat, lon) in coords.items():
            if code not in counts:
                ax.scatter(lon, lat, c="k", marker="x", alpha=0.5, s=75,
                           zorder=5)

        _set_map_axes(ax, extent)
        ax.set_title(f"{len(codes)} {key} stations\n"
                     f"N_measurements=[{min(cnts)}, {max(cnts)}]")
        f.savefig(os.path.join(output_dir, f"{key}_counts.png"))
        plt.close(f)


def _plot_source_coverage(ax, src_code, kernels, coords, extent=None):
    """
    Plot one source station and the receivers that have data for it

    :type ax: matplotlib.axes.Axes
    :param ax: axis to plot on
    :type src_code: str
    :param src_code: 'NET_STA' code of the source station
    :type kernels: dict
    :param kernels: {'ZZ' or 'TT': list of 'NET_STA' receiver codes}
    :type coords: dict
    :param coords: {'NET_STA': (lat, lon)} from :func:`read_station_coords`
    :type extent: list of float or None
    :param extent: [lon_min, lon_max, lat_min, lat_max], None for automatic
    :rtype: int
    :return: total number of measurements plotted
    """
    src_lat, src_lon = coords[src_code]
    styles = {"ZZ": dict(c="b", marker="v"), "TT": dict(c="g", marker="^")}

    for lat, lon in coords.values():
        ax.scatter(lon, lat, c="k", marker="s", zorder=3, alpha=0.25)
    ax.scatter(src_lon, src_lat, c="r", marker="o", zorder=6)

    n = 0
    for kernel, rcv_codes in sorted(kernels.items()):
        rcv_codes = [code for code in rcv_codes if code in coords]
        if not rcv_codes:
            continue
        lats, lons = zip(*[coords[code] for code in rcv_codes])
        ax.scatter(lons, lats, zorder=5, label=f"{kernel} ({len(lats)})",
                   **styles.get(kernel, {}))
        for lat, lon in zip(lats, lons):
            ax.plot([src_lon, lon], [src_lat, lat], c="k", ls="-",
                    alpha=0.25, zorder=4)
        n += len(lats)

    _set_map_axes(ax, extent)
    ax.legend(loc="upper right")
    ax.set_title(f"{src_code.replace('_', '.')}; N={n}")
    return n


def plot_measurement_coverage(data_dir, coords, output_dir, extent=None):
    """
    Plot source-receiver paths for each source, receiver and doubled source

    Writes to subdirectories 'sources', 'receivers' and 'doubles' of
    ``output_dir``. The 'doubles' figures show stations sharing a code
    across networks side-by-side, to visually decide which one to keep.

    :type data_dir: str
    :param data_dir: output directory of :func:`prep`
    :type coords: dict
    :param coords: {'NET_STA': (lat, lon)} from :func:`read_station_coords`
    :type output_dir: str
    :param output_dir: directory to save figures to
    :type extent: list of float or None
    :param extent: [lon_min, lon_max, lat_min, lat_max], None for automatic
    """
    src_paths = defaultdict(lambda: defaultdict(list))
    rcv_paths = defaultdict(list)
    for src_code, kernel, rcv_code in iter_prepped_files(data_dir):
        src_paths[src_code][kernel].append(rcv_code)
        rcv_paths[rcv_code].append(src_code)

    missing = (set(src_paths) | set(rcv_paths)) - set(coords)
    if missing:
        print(f"WARNING: {len(missing)} stations with data have no "
              f"coordinates and are not plotted: {sorted(missing)}")

    for subdir in ["sources", "receivers", "doubles"]:
        os.makedirs(os.path.join(output_dir, subdir), exist_ok=True)

    for src_code, kernels in src_paths.items():
        if src_code not in coords:
            continue
        f, ax = plt.subplots()
        n = _plot_source_coverage(ax, src_code, kernels, coords, extent)
        f.savefig(os.path.join(output_dir, "sources",
                               f"{n}_{src_code}.png"))
        plt.close(f)

    for rcv_code, src_codes in rcv_paths.items():
        if rcv_code not in coords:
            continue
        rcv_lat, rcv_lon = coords[rcv_code]
        src_codes = [code for code in src_codes if code in coords]
        f, ax = plt.subplots()
        ax.scatter(rcv_lon, rcv_lat, c="r", marker="v", zorder=6)
        for src_code in src_codes:
            src_lat, src_lon = coords[src_code]
            ax.scatter(src_lon, src_lat, c="g", marker="o", zorder=5)
            ax.plot([src_lon, rcv_lon], [src_lat, rcv_lat], c="k", ls="-",
                    alpha=0.25, zorder=4)
        _set_map_axes(ax, extent)
        n = len(src_codes)
        ax.set_title(f"{rcv_code.replace('_', '.')}; N={n}")
        f.savefig(os.path.join(output_dir, "receivers",
                               f"{n}_{rcv_code}_rcv.png"))
        plt.close(f)

    src_count = {code: sum(len(v) for v in kernels.values())
                 for code, kernels in src_paths.items() if code in coords}
    for sta, codes in find_doubled_sources(src_count).items():
        f, axs = plt.subplots(1, len(codes), figsize=(6.4 * len(codes), 4.8))
        for ax, (src_code, _) in zip(axs, codes):
            _plot_source_coverage(ax, src_code, src_paths[src_code], coords,
                                  extent)
        f.savefig(os.path.join(output_dir, "doubles", f"{sta}.png"))
        plt.close(f)


# =============================================================================
#                                COMMAND LINE
# =============================================================================
def parse_args(argv=None):
    """
    Parse command line arguments

    :type argv: list of str or None
    :param argv: arguments to parse, None to use sys.argv
    :rtype: argparse.Namespace
    :return: parsed arguments
    """
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("inspect", help="print SAC headers of input files")
    p.add_argument("input", help="extracted Zenodo directory")
    p.add_argument("-n", "--nfiles", type=int, default=3,
                   help="files to print per stack directory")
    p.add_argument("--old-origin", default=OLD_ORIGIN,
                   help="zero-lag time of the original data")

    p = sub.add_parser("prep", help="write a modified copy of the dataset")
    p.add_argument("input", help="extracted Zenodo directory")
    p.add_argument("output", help="new directory for the modified data")
    p.add_argument("--type", default="hyp", choices=["hyp", "ell"],
                   help="data type to keep")
    p.add_argument("--sources", default=None,
                   help="STATIONS file of source stations to keep")
    p.add_argument("--stations", default=None,
                   help="STATIONS file of receiver stations to keep")
    p.add_argument("--old-origin", default=OLD_ORIGIN,
                   help="zero-lag time of the original data")
    p.add_argument("--new-origin", default=NEW_ORIGIN,
                   help="zero-lag time of the output data")
    p.add_argument("--length", type=float, default=LENGTH_S,
                   help="output trace length in seconds")
    p.add_argument("--overwrite", action="store_true",
                   help="re-process files that already exist in output")
    p.add_argument("--dry-run", action="store_true",
                   help="count files that would be written, write nothing")

    p = sub.add_parser("select", help="write thresholded STATIONS files")
    p.add_argument("data", help="output directory of 'prep'")
    p.add_argument("stations_all", help="STATIONS file of all stations")
    p.add_argument("--src-threshold", type=int, default=86,
                   help="minimum measurements per source")
    p.add_argument("--rcv-threshold", type=int, default=5,
                   help="minimum measurements per receiver")
    p.add_argument("-o", "--output", default=".",
                   help="directory to write files to")

    p = sub.add_parser("plot", help="make station coverage figures")
    p.add_argument("data", help="output directory of 'prep'")
    p.add_argument("stations_all", help="STATIONS file of all stations")
    p.add_argument("--src-threshold", type=int, default=86,
                   help="minimum measurements per source for count maps")
    p.add_argument("--rcv-threshold", type=int, default=5,
                   help="minimum measurements per receiver for count maps")
    p.add_argument("--extent", type=float, nargs=4, default=None,
                   metavar=("LONMIN", "LONMAX", "LATMIN", "LATMAX"),
                   help="map extent, e.g., -168 -140 64.5 72")
    p.add_argument("-o", "--output", default="coverage_figures",
                   help="directory to save figures to")

    return parser.parse_args(argv)


def main(argv=None):
    """
    Run the subcommand chosen on the command line

    :type argv: list of str or None
    :param argv: arguments to parse, None to use sys.argv
    """
    args = parse_args(argv)
    if args.command == "inspect":
        inspect(args.input, args.nfiles, args.old_origin)
    elif args.command == "prep":
        prep(args.input, args.output, data_type=args.type,
             sources_file=args.sources, stations_file=args.stations,
             old_origin=args.old_origin, new_origin=args.new_origin,
             length_s=args.length, overwrite=args.overwrite,
             dry_run=args.dry_run)
    elif args.command == "select":
        select_stations(args.data, args.stations_all, args.src_threshold,
                        args.rcv_threshold, args.output)
    elif args.command == "plot":
        os.makedirs(args.output, exist_ok=True)
        coords = read_station_coords(args.stations_all)
        src_count, rcv_count = count_source_receiver_hits(
            args.data, args.src_threshold, args.rcv_threshold)
        plot_source_receiver_hits(src_count, rcv_count, coords, args.output,
                                  args.extent)
        plot_measurement_coverage(args.data, coords, args.output,
                                  args.extent)


if __name__ == "__main__":
    sys.exit(main())
