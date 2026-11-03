#!/usr/bin/env python
"""
Batch driver for psd_tools.py that replaces run_psd_viz.sh.

Scans the MiniSEED directory (DATA) for every available channel and runs
three steps:

1. ``compute``: build one combined PPSD NPZ per channel in SPECTRA
2. ``pdf``: plot a PPSD PDF figure for each NPZ into FIGURES
3. ``compare``: plot per-station comparisons, grouping each main station
   (e.g., AK.MDM) with its sub-stations (e.g., XX.MDMC, XX.MDMU, XX.MDMW),
   one figure per component

Every task's outcome is recorded in a JSON state file. On re-runs, tasks that
succeeded (and whose output still exists) are skipped, and tasks that failed
are also skipped so a broken channel is not re-attempted every time. Use
``--retry-failed`` to try failed tasks again, or ``--force`` to redo
everything.

Example::

    python run_psd_viz.py                     # run all steps
    python run_psd_viz.py --steps compute     # only make spectra
    python run_psd_viz.py --stations MDM,H23K # restrict to station groups
    python run_psd_viz.py --dry-run           # show what would run
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path


# Default paths, mirroring the original shell script
PYTHON = "/Users/prof/micromamba/envs/psd_viz/bin/python"
SCRIPT = "/Users/prof/repos/psd_viz/psd_tools.py"
PPSDS = "SPECTRA"
DATA = "DATA"
FIGS = "FIGURES"
LOGS = "logs"
XX_XML = "ISHTAR.xml"
STATE_FILE = "run_psd_viz_state.json"

# Expected file names: NET.STA.LOC.CHA.YEAR.JDAY, e.g. AK.MDM..BHZ.2025.190
MSEED_RE = re.compile(
    r"^(?P<net>[A-Z0-9]+)\.(?P<sta>[A-Z0-9]+)\.(?P<loc>[A-Z0-9]*)\."
    r"(?P<cha>[A-Z0-9]{3})\.(?P<year>\d{4})\.(?P<jday>\d{3})$"
)

# Only process seismic channels (instrument code 'H'); this excludes
# non-seismic channels such as LWD/LWS wind data
DEFAULT_CHANNEL_RE = r"^.H[ZNE12]$"

# Network that holds the sub-stations of each main station
SUB_NETWORK = "XX"

STEPS = ("compute", "pdf", "compare")


@dataclass
class Channel:
    """
    A single SNCL found in the data directory, with its data time span.

    :type net: str
    :param net: network code
    :type sta: str
    :param sta: station code
    :type loc: str
    :param loc: location code, may be empty
    :type cha: str
    :param cha: channel code
    :type first_day: datetime
    :param first_day: UTC midnight of the first day with data
    :type last_day: datetime
    :param last_day: UTC midnight of the last day with data
    :type nfiles: int
    :param nfiles: number of day files found for this channel
    """
    net: str
    sta: str
    loc: str
    cha: str
    first_day: datetime
    last_day: datetime
    nfiles: int

    @property
    def sncl(self):
        """Dot-separated SNCL code, e.g. AK.MDM..BHZ"""
        return f"{self.net}.{self.sta}.{self.loc}.{self.cha}"

    @property
    def npz_name(self):
        """File name psd_tools.py gives the combined NPZ for this channel"""
        parts = [self.net, self.sta]
        if self.loc:
            parts.append(self.loc)
        parts.append(self.cha)
        return "_".join(parts) + "_all.npz"

    @property
    def pdf_name(self):
        """File name psd_tools.py gives the PDF figure for this channel"""
        return self.npz_name.replace("_all.npz", "_pdf.png")

    @property
    def component(self):
        """Single-letter component code, e.g. Z"""
        return self.cha[-1]


class State:
    """
    Persistent record of task outcomes so finished or failed work is not
    repeated on later runs.

    :type path: pathlib.Path
    :param path: JSON file that stores the state
    """
    def __init__(self, path):
        self.path = Path(path)
        self.tasks = {}
        if self.path.exists():
            try:
                self.tasks = json.loads(self.path.read_text())
            except json.JSONDecodeError:
                print(f"WARNING: could not parse {self.path}, starting fresh")

    def get(self, key):
        """
        Return the stored record for a task, or None.

        :type key: str
        :param key: unique task identifier
        :rtype: dict or None
        :return: stored task record
        """
        return self.tasks.get(key)

    def set(self, key, status, **info):
        """
        Record a task outcome and write the state file immediately so that an
        interrupted run keeps everything finished so far.

        :type key: str
        :param key: unique task identifier
        :type status: str
        :param status: 'ok' or 'failed'
        :type info: dict
        :param info: extra fields to store, e.g. output path or error text
        """
        self.tasks[key] = {
            "status": status,
            "time": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            **info,
        }
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps(self.tasks, indent=2, sort_keys=True))
        tmp.replace(self.path)


class Runner:
    """
    Runs psd_tools.py subcommands with skip/fail bookkeeping.

    :type args: argparse.Namespace
    :param args: parsed command line arguments
    """
    def __init__(self, args):
        self.args = args
        self.state = State(Path(args.logs) / STATE_FILE)
        self.counts = defaultdict(int)

    def should_skip(self, key, outputs, members=None):
        """
        Decide whether a task can be skipped based on its prior outcome.

        :type key: str
        :param key: unique task identifier
        :type outputs: list of pathlib.Path
        :param outputs: files the task is expected to produce
        :type members: list of str
        :param members: optional input list; a change forces a re-run
        :rtype: bool
        :return: True if the task should not be run
        """
        if self.args.force:
            return False
        record = self.state.get(key)
        if record is None:
            # No record, but outputs from earlier (pre-state) runs count
            if members is None and outputs and all(
                    p.exists() for p in outputs):
                print(f"  skip  {key}: output exists")
                self.counts["skipped"] += 1
                return True
            return False
        if members is not None and record.get("members") != members:
            print(f"  redo  {key}: inputs changed")
            return False
        if record["status"] == "ok":
            if all(p.exists() for p in outputs):
                print(f"  skip  {key}: already done")
                self.counts["skipped"] += 1
                return True
            print(f"  redo  {key}: output missing")
            return False
        if record["status"] == "failed" and not self.args.retry_failed:
            print(f"  skip  {key}: failed previously "
                  f"({record.get('error', 'unknown error')}); "
                  f"use --retry-failed to try again")
            self.counts["skipped_failed"] += 1
            return True
        return False

    def run(self, key, cmd, outputs, log_name, members=None):
        """
        Run one psd_tools.py command, then verify its outputs and record the
        result in the state file.

        :type key: str
        :param key: unique task identifier
        :type cmd: list of str
        :param cmd: command to execute
        :type outputs: list of pathlib.Path
        :param outputs: files that must exist afterwards for success
        :type log_name: str
        :param log_name: base name of the log file for this task
        :type members: list of str
        :param members: optional input list stored with the record
        :rtype: bool
        :return: True if the task succeeded
        """
        print(f"  run   {key}")
        if self.args.dry_run:
            print("        " + " ".join(cmd))
            return True

        log_path = Path(self.args.logs) / "run_psd_viz" / f"{log_name}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        extra = {} if members is None else {"members": members}
        try:
            with open(log_path, "w") as log:
                log.write(" ".join(cmd) + "\n\n")
                log.flush()
                proc = subprocess.run(cmd, stdout=log,
                                      stderr=subprocess.STDOUT)
        except (OSError, KeyboardInterrupt) as e:
            if isinstance(e, KeyboardInterrupt):
                # Do not mark as failed; it should run again next time
                print("\nInterrupted; this task will be re-run next time")
                raise
            error = f"could not execute: {e}"
            proc = None

        if proc is not None:
            missing = [str(p) for p in outputs if not p.exists()]
            if proc.returncode != 0:
                error = f"exit code {proc.returncode}"
            elif missing:
                error = f"no output produced: {', '.join(missing)}"
            else:
                self.state.set(key, "ok",
                               outputs=[str(p) for p in outputs], **extra)
                self.counts["ok"] += 1
                return True

        print(f"  FAIL  {key}: {error} (see {log_path})")
        self.state.set(key, "failed", error=error, log=str(log_path),
                       **extra)
        self.counts["failed"] += 1
        return False


def scan_data(data_dir, channel_re):
    """
    Find every channel in the data directory and the days it covers.

    :type data_dir: str
    :param data_dir: directory of day-long MiniSEED files
    :type channel_re: str
    :param channel_re: regular expression channel codes must match
    :rtype: list of Channel
    :return: channels sorted by network, station, channel
    """
    pattern = re.compile(channel_re)
    days = defaultdict(list)
    for name in os.listdir(data_dir):
        match = MSEED_RE.match(name)
        if not match or not pattern.match(match["cha"]):
            continue
        key = (match["net"], match["sta"], match["loc"], match["cha"])
        day = datetime(int(match["year"]), 1, 1, tzinfo=timezone.utc)
        days[key].append(day + timedelta(days=int(match["jday"]) - 1))

    channels = [
        Channel(net, sta, loc, cha, min(d), max(d), len(d))
        for (net, sta, loc, cha), d in days.items()
    ]
    return sorted(channels, key=lambda c: (c.net, c.sta, c.cha))


def main_station(channel, main_names):
    """
    Map a channel to the main station it belongs to. Sub-stations in the
    XX network are named after their main station plus one suffix letter,
    e.g. MDMC, MDMU, MDMW all belong to MDM.

    :type channel: Channel
    :param channel: channel to classify
    :type main_names: set of str
    :param main_names: station codes of all main (non-XX) stations
    :rtype: str
    :return: main station code
    """
    if channel.net != SUB_NETWORK:
        return channel.sta
    # Prefer the longest main station code that prefixes this station
    for name in sorted(main_names, key=len, reverse=True):
        if channel.sta.startswith(name) and channel.sta != name:
            return name
    # No main station in the data; fall back to dropping the suffix letter
    return channel.sta[:-1]


def stationxml_for(channel, root):
    """
    Return the StationXML file holding the response for a channel.

    :type channel: Channel
    :param channel: channel needing a response
    :type root: pathlib.Path
    :param root: directory that holds the StationXML files
    :rtype: pathlib.Path
    :return: path to the StationXML file
    """
    if channel.net == SUB_NETWORK:
        return root / XX_XML
    return root / f"{channel.net}.{channel.sta}.xml"


def step_compute(runner, channels):
    """
    Create one combined PPSD NPZ per channel.

    The NPZ is written into a scratch directory and only moved into place
    once psd_tools.py finishes, so a crash never leaves a partial NPZ that
    would be mistaken for a finished one.

    :type runner: Runner
    :param runner: task runner
    :type channels: list of Channel
    :param channels: channels to process
    """
    args = runner.args
    spectra = Path(args.spectra)
    root = Path(args.root)
    print(f"\n=== Step 1: compute PPSDs ({len(channels)} channels) ===")
    for chan in channels:
        key = f"compute:{chan.sncl}"
        final = spectra / chan.npz_name
        if runner.should_skip(key, [final]):
            continue

        xml = stationxml_for(chan, root)
        if not xml.exists():
            print(f"  FAIL  {key}: missing StationXML {xml}")
            if not args.dry_run:
                runner.state.set(key, "failed",
                                 error=f"missing StationXML {xml}")
                runner.counts["failed"] += 1
            continue

        tmpdir = spectra / f".tmp_{chan.npz_name[:-4]}"
        start = chan.first_day.strftime("%Y-%m-%d")
        end = (chan.last_day + timedelta(days=1)).strftime("%Y-%m-%d")
        cmd = [
            args.python, args.script, "compute",
            start, end, chan.net, chan.sta, chan.cha,
            "--loc", chan.loc or "*",
            "--outdir", str(tmpdir),
            "--mseed-dir", args.data,
            "--mseed-pattern", f"{chan.sncl}.*",
            "--stationxml", str(xml),
            "--num-proc", str(args.num_proc),
        ]
        if not args.dry_run:
            shutil.rmtree(tmpdir, ignore_errors=True)
            tmpdir.mkdir(parents=True)
        ok = runner.run(key, cmd, [] if args.dry_run else
                        [tmpdir / chan.npz_name],
                        log_name=f"compute_{chan.npz_name[:-8]}")
        if args.dry_run:
            continue
        if ok:
            (tmpdir / chan.npz_name).replace(final)
            # Point the stored output at the final location
            runner.state.set(key, "ok", outputs=[str(final)])
        shutil.rmtree(tmpdir, ignore_errors=True)


def step_pdf(runner, channels):
    """
    Plot a PPSD PDF figure for each channel that has an NPZ.

    :type runner: Runner
    :param runner: task runner
    :type channels: list of Channel
    :param channels: channels to plot
    """
    args = runner.args
    spectra = Path(args.spectra)
    figs = Path(args.figures)
    print(f"\n=== Step 2: plot PPSD PDFs ({len(channels)} channels) ===")
    for chan in channels:
        key = f"pdf:{chan.sncl}"
        npz = spectra / chan.npz_name
        if not npz.exists():
            print(f"  skip  {key}: no NPZ (compute step missing or failed)")
            runner.counts["no_input"] += 1
            continue
        out = figs / chan.pdf_name
        if runner.should_skip(key, [out]):
            continue
        cmd = [
            args.python, args.script, "pdf", str(npz),
            "--cmap", "cividis",
            "--color-by", "probability",
            "--show-percentiles",
            "--domain", "frequency",
            "--outdir", str(figs),
        ]
        runner.run(key, cmd, [out], log_name=f"pdf_{chan.npz_name[:-8]}")


def step_compare(runner, channels):
    """
    Plot comparison figures per main station and component, combining the
    main station with all of its XX sub-stations.

    :type runner: Runner
    :param runner: task runner
    :type channels: list of Channel
    :param channels: channels to group and compare
    """
    args = runner.args
    spectra = Path(args.spectra).resolve()
    figs = Path(args.figures)
    main_names = {c.sta for c in channels if c.net != SUB_NETWORK}

    groups = defaultdict(list)
    for chan in channels:
        groups[(main_station(chan, main_names), chan.component)].append(chan)

    print(f"\n=== Step 3: station comparisons ({len(groups)} groups) ===")
    for (main, comp), members in sorted(groups.items()):
        key = f"compare:{main}.{comp}"
        npzs = sorted(spectra / c.npz_name for c in members
                      if (spectra / c.npz_name).exists())
        if len(npzs) < 2:
            print(f"  skip  {key}: fewer than 2 NPZs available")
            runner.counts["no_input"] += 1
            continue
        names = [p.name for p in npzs]
        out = figs / f"{main}_{comp}_compare.png"
        if runner.should_skip(key, [out], members=names):
            continue

        # psd_tools.py compares every NPZ in a directory, so gather this
        # group's files into their own directory of symlinks
        group_dir = spectra / f"{main}_{comp}_compare"
        if not args.dry_run:
            shutil.rmtree(group_dir, ignore_errors=True)
            group_dir.mkdir(parents=True)
            for npz in npzs:
                (group_dir / npz.name).symlink_to(
                    os.path.relpath(npz, group_dir))
        cmd = [
            args.python, args.script, "compare", str(group_dir),
            "--outfile", str(out),
            "--shade-percentiles",
            "--domain", "frequency",
            "--title", f"{main} {comp}-component comparison",
        ]
        runner.run(key, cmd, [out], log_name=f"compare_{main}_{comp}",
                   members=names)


def parse_args():
    """
    Parse command line arguments.

    :rtype: argparse.Namespace
    :return: parsed arguments
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--steps", default=",".join(STEPS),
                        help="Comma list of steps to run: "
                             f"{','.join(STEPS)}. Default: all")
    parser.add_argument("--stations", default=None,
                        help="Comma list of main stations to process, "
                             "e.g. MDM,H23K (includes sub-stations)")
    parser.add_argument("--channel-regex", default=DEFAULT_CHANNEL_RE,
                        help="Regex channel codes must match. "
                             f"Default: {DEFAULT_CHANNEL_RE}")
    parser.add_argument("--num-proc", type=int, default=1,
                        help="Worker processes for psd_tools.py compute")
    parser.add_argument("--retry-failed", action="store_true",
                        help="Re-attempt tasks that failed before")
    parser.add_argument("--force", action="store_true",
                        help="Re-run every task, even if already done")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would run without running it")
    parser.add_argument("--python", default=PYTHON,
                        help=f"Python interpreter. Default: {PYTHON}")
    parser.add_argument("--script", default=SCRIPT,
                        help=f"psd_tools.py path. Default: {SCRIPT}")
    parser.add_argument("--root", default=".",
                        help="Directory holding the StationXML files")
    parser.add_argument("--data", default=DATA, help="MiniSEED directory")
    parser.add_argument("--spectra", default=PPSDS,
                        help="NPZ output directory")
    parser.add_argument("--figures", default=FIGS,
                        help="Figure output directory")
    parser.add_argument("--logs", default=LOGS,
                        help="Directory for logs and the state file")
    args = parser.parse_args()

    args.steps = [s.strip() for s in args.steps.split(",") if s.strip()]
    bad = set(args.steps) - set(STEPS)
    if bad:
        parser.error(f"unknown step(s): {', '.join(sorted(bad))}")
    if args.stations:
        args.stations = {s.strip() for s in args.stations.split(",")}
    return args


def main():
    """Scan the data directory and run the requested steps"""
    args = parse_args()
    for path in (args.python, args.script, args.data):
        if not Path(path).exists():
            sys.exit(f"ERROR: {path} does not exist")
    for path in (args.spectra, args.figures, args.logs):
        Path(path).mkdir(parents=True, exist_ok=True)

    channels = scan_data(args.data, args.channel_regex)
    if args.stations:
        main_names = {c.sta for c in channels if c.net != SUB_NETWORK}
        channels = [c for c in channels
                    if main_station(c, main_names) in args.stations]
    if not channels:
        sys.exit("ERROR: no matching channels found")

    print(f"Found {len(channels)} channels in {args.data}:")
    for chan in channels:
        print(f"  {chan.sncl:<18} {chan.first_day:%Y-%m-%d} to "
              f"{chan.last_day:%Y-%m-%d} ({chan.nfiles} files)")

    runner = Runner(args)
    try:
        if "compute" in args.steps:
            step_compute(runner, channels)
        if "pdf" in args.steps:
            step_pdf(runner, channels)
        if "compare" in args.steps:
            step_compare(runner, channels)
    except KeyboardInterrupt:
        sys.exit(130)

    print("\n=== Summary ===")
    for name in ("ok", "skipped", "failed", "skipped_failed", "no_input"):
        print(f"  {name:<15} {runner.counts[name]}")
    print(f"State file: {runner.state.path}")
    if runner.counts["failed"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
