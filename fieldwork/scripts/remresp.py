"""
Remove response from Magseis-Fairfield ZLand 1C and 3C nodal data, or 
SmartSolo IGU-BD3C-5 nodes. Builds responses directly from Nominal Response
Library (NRL) URLs based on user-defined setup parameters. Returns new MSEED 
files with the same name containing response removed data.

You must choose which instrument your data is from with the <choice> parameter, 
this will help select the available choices for parameters. See below.

.. rubric::
        
    $ python remresp.py <choice> <files> <flags> 

    where choices are currently 'fairfield' and 'smartsolo' e.g.,

    $ python remresp.py fairfield XX.101..DHZ.2026.* \
        --sample_rate 250 --output_units count --pre_amp_gain 18 \
        --filter_phase LP --dc_filter Off --pre_filt .001 .005 120 125 \
        --output VEL --water_level 60

    # Turn off default options including tapering, demean, and writing inv
    $ python remresp.py smartsolo XX.101..DHZ.2026.101 --taper_fraction 0 \
            --write_inv '' --zero_mean_off

.. rubric:: Offline workflow with a comprehensive StationXML

    Rather than removing response directly, `--make_inv` builds a single
    StationXML file containing every channel found in the input files, with
    coordinates taken from a user-defined text/CSV file (see `read_coords`
    for the format). Response removal is skipped in this mode.

    $ python remresp.py fairfield XX.*.DHZ.2026.* -r 250 -p 18 -d Off \
        --make_inv --coords stations.csv --write_inv XX_inventory.xml

.. note:: Mixed sample rates

    If `--sample_rate` is not given, the sample rate of each file is read
    from its header and the matching NRL response is fetched (once per unique
    sample rate). This allows datasets with multiple sample rates to be
    processed together. With `--make_inv`, a station that changed sample
    rate receives one channel epoch per sample rate.

    The resulting StationXML can then be used to remove response without
    accessing NRL webservices using the 'stationxml' choice

    $ python remresp.py stationxml XX.*.DHZ.2026.* --inv XX_inventory.xml \
        --pre_filt .001 .005 120 125 --output VEL

.. note:: count or mV kludge

    The warning message below is given if our input data is in units 'mV' 
    because this does not match ObsPy's internal mapping dictionaries. This is
    fine because we only want the gain that is added from this stage so we 
    ignore the warning message to avoid confusion.

    UserWarning: The unit 'MV' is not known to ObsPy. It will be passed in to 
    evalresp as 'undefined'. This should result in evalresp using the response 
    as is, without adding any integration or differentiation and the 'output'
    parameter (here: 'VEL') not having any effect.`

.. note:: Fairfield Sign Convention

    If Fairfield nodal data were converted using the `fcnt2mseed.py` script, 
    then metadata should define `dip=-90` to maintain +Z up orientation that 
    is enforced by the ObsPy read function. This is done by default here.
"""
import os
import re
import argparse
import warnings
from copy import deepcopy
from obspy import read, read_inventory, UTCDateTime
from obspy.core.inventory import Inventory, Network, Station, Channel

# Available parameters choices for each node type
PARAMETERS = {
    # Magseis Fairfield ZLand 3C
    # https://ds.iris.edu/ds/nrl/datalogger/magseisfairfield/zlandgen2/
    # https://ds.iris.edu/ds/nrl/sensor/magseisfairfield/zlandgen2sensor/
    "fairfield": {
        "sample_rate": (250, 500, 1000, 2000),
        "output_units": ("count", "mV"),
        "pre_amp_gain": (0, 6, 12, 18, 24, 30, 36),
        "filter_phase": ("LP", "MP"),
        "dc_filter": ("1", "Off"),
        },
    # SmartSolo IGU-BD3C-5
    # https://ds.iris.edu/ds/nrl/sensor/dtcc/dt-solo-bb/
    # https://ds.iris.edu/ds/nrl/datalogger/dtcc/smartsolo-igu-bd3c-5/
    "smartsolo": {
        "sample_rate": (50, 100, 125, 250, 500, 1000, 2000, 4000),
        "output_units": ("count", "mV"),
        "pre_amp_gain": (0, 6),
        "filter_phase": ("LP", "MP"),
        "dc_filter": ("1", "DC", "Off")
        }
    }
BASE_NRL_URL = "https://service.earthscope.org/irisws/nrl/1/combine?"

# Format code -> (attribute name, regex). Code fields have flexible widths so
# 'S' or 'SSS' both match 3--5 char station names; date fields have fixed
# widths so that delimiter-free dates (e.g., YYYYMMDD) can still be split
FIELDS = {
    "N": ("network", r"[A-Za-z0-9]{1,2}"),
    "S": ("station", r"[A-Za-z0-9]{1,5}"),
    "L": ("location", r"[A-Za-z0-9]{0,2}"),
    "C": ("channel", r"[A-Za-z0-9]{3}"),
    "Y": ("year", r"\d{4}"),
    "M": ("month", r"\d{2}"),
    "D": ("day", r"\d{2}"),
    "J": ("julday", r"\d{3}"),
}

# Accepted header names (lowercase) for each column of the coordinate file
COORD_COLUMNS = {
    "network": ("network", "net"),
    "station": ("station", "sta", "name"),
    "latitude": ("latitude", "lat"),
    "longitude": ("longitude", "lon", "long"),
    "elevation": ("elevation", "elev", "elv"),
    "depth": ("depth",),
}

# Channel orientation code -> (azimuth, dip) following SEED convention, where
# dip=-90 is +Z up. See the Fairfield sign convention note above
ORIENTATIONS = {
    "Z": (0.0, -90.0),
    "N": (0.0, 0.0),
    "1": (0.0, 0.0),
    "E": (90.0, 0.0),
    "2": (90.0, 0.0),
}


def parse_args():
    """
    Parse command line arguments. The first positional argument selects the
    instrument, which sets the valid choices for its response parameters.

    :rtype: argparse.Namespace
    :return: parsed arguments; `choice` holds the instrument name
    """
    # Arguments shared by every instrument
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("fids", nargs="+", help="required, file ID(s)")

    # Arguments fed into `trace.remove_response()`
    common.add_argument("-f", "--pre_filt", nargs=4, type=float,
                        default=None, metavar=("F1", "F2", "F3", "F4"),
                        help="optional pre filter corners [Hz]")
    common.add_argument("--water_level", default=60, type=float,
                        help="apply a water level in the response removal")
    common.add_argument("--zero_mean_off", action="store_true",
                        help="by default response removal zeroes the mean, you "
                             "can turn it off with this flag (i.e., you will "
                             "not demean data")
    common.add_argument("--taper_fraction", type=float, default=0.05,
                        help="apply a taper during response removal. If 0, "
                             "no taper is applied")
    common.add_argument("--output", type=str, default="VEL", 
                        choices=("DISP", "VEL", "ACC"),
                        help="output ground motion quantity")

    # File saving
    common.add_argument("--save", default="./response_removed",
                        help="where to save the newly created files")
    common.add_argument("--overwrite", action="store_true",
                        help="overwrite any existing files in `save`")
    common.add_argument("--rename", nargs="?", default=None,
                        const="NN.SSS.LL.CCC.YYYY.JJJ", metavar="FMT",
                        help="set internal trace stats (network, station, "
                             "location, channel) from the file name, for "
                             "when raw stats do not match user-defined "
                             "naming. Optionally provide a format string "
                             "FMT describing the file name using codes "
                             "N=network, S=station, L=location, C=channel, "
                             "Y=year, M=month, D=day, J=julian day. Other "
                             "characters are treated as delimiters and file "
                             "extensions are ignored. If no FMT is given, "
                             "defaults to 'NN.SSS.LL.CCC.YYYY.JJJ'. Code "
                             "letter counts are not enforced (S matches 1-5 "
                             "chars) but date codes are fixed width. Wrap "
                             "literal text containing code letters in square "
                             "brackets. e.g., --rename YYYYMMDD.NN.SSS.LL.CCC "
                             "or --rename 'NN.SSS.LL.CCC.YYYY.JJJ[_SLICED]'")
    common.add_argument("--write_inv", type=str, default="inventory.xml",
                        help="filename to write the response out as "
                             "a StationXML file for archival. This will write "
                             "the inventory for the first data file only, "
                             "or for all files if `--make_inv` is set. "
                             "Set as empty string to toggle off")

    # Instrument specific parameters
    parser = argparse.ArgumentParser(
        description="Remove nominal NRL response from node data")
    sub = parser.add_subparsers(dest="choice", required=True,
                                metavar="{fairfield,smartsolo,stationxml}")

    # Remove response with an existing StationXML file, no NRL access needed
    p = sub.add_parser("stationxml", parents=[common],
                       help="use an existing StationXML file, e.g., one "
                            "created with `--make_inv`")
    p.add_argument("-i", "--inv", type=str, required=True,
                   help="path to StationXML file containing responses for "
                        "all channels in `fids`")

    helps = {"fairfield": "Magseis Fairfield ZLand 3C",
             "smartsolo": "SmartSolo IGU-BD3C-5"}
    for name, opts in PARAMETERS.items():
        p = sub.add_parser(name, parents=[common], help=helps[name])
        p.add_argument("-r", "--sample_rate", type=int, default=None,
                       choices=opts["sample_rate"],
                       help="final sample rate [Hz]. If not given, the "
                            "sample rate is read from each file header")
        p.add_argument("-u", "--output_units", default="mV",
                       choices=opts["output_units"],
                       help="units of the data on disk")
        p.add_argument("-p", "--pre_amp_gain", type=int, required=True,
                       choices=opts["pre_amp_gain"],
                       help="preamp gain [dB]")
        p.add_argument("-t", "--filter_phase", default="LP",
                       choices=opts["filter_phase"],
                       help="LP=linear phase, MP=minimum phase")
        p.add_argument("-d", "--dc_filter", required=True,
                       choices=opts["dc_filter"],
                       help="low-cut filter setting")
        p.add_argument("--make_inv", action="store_true",
                       help="instead of removing response, write a single "
                            "StationXML file (to `write_inv`) containing "
                            "all channels in `fids`. Requires `--coords`")
        p.add_argument("-c", "--coords", type=str, default=None,
                       help="text or CSV file with station coordinates, "
                            "used with `--make_inv`. Columns are "
                            "'station latitude longitude [elevation]' or "
                            "set by a header row, see `read_coords`")

    args = parser.parse_args()
    if getattr(args, "make_inv", False):
        if not args.coords:
            parser.error("`--make_inv` requires `--coords`")
        if not args.write_inv:
            parser.error("`--make_inv` requires a `--write_inv` filename")

    return args


def parse_fname(fname, fmt="NN.SSS.LL.CCC.YYYY.JJJ"):
    """
    Split a file name into its components based on a user-defined format
    string. Runs of a format code (e.g., 'SSS') count as one field, all other
    characters are treated as literal delimiters. Text in square brackets is
    matched literally, even if it contains format codes (e.g., '[_SLICED]').
    Repeated codes are matched but only the first occurrence is kept. Any
    trailing file extension (e.g., .ms, .mseed, .MSEED) is ignored.

    :type fname: str
    :param fname: file name or path to parse, only the basename is used
    :type fmt: str
    :param fmt: format string built from codes N=network, S=station,
        L=location, C=channel, Y=year, M=month, D=day, J=julian day
    :rtype: dict
    :return: parsed components keyed by attribute name, e.g., 'station'
    :raises ValueError: if the file name does not match the format
    """
    pattern = ""
    seen = set()
    # Tokens are either [literal text] or a run of one repeated character
    for token in re.finditer(r"\[(?P<literal>[^\]]*)\]|(?P<run>(.)\3*)", fmt):
        if token.group("literal") is not None:
            pattern += re.escape(token.group("literal"))
            continue
        run = token.group("run")
        if run[0] in FIELDS:
            name, regex = FIELDS[run[0]]
            # Only capture the first occurrence of a code; repeats are matched
            # but discarded so that group names are not redefined
            if name in seen:
                pattern += f"(?:{regex})"
            else:
                pattern += f"(?P<{name}>{regex})"
                seen.add(name)
        else:
            pattern += re.escape(run)
    pattern += r"(?:\.[A-Za-z0-9]+)?"  # optional file extension

    match = re.fullmatch(pattern, os.path.basename(fname))
    if match is None:
        raise ValueError(f"'{fname}' does not match format '{fmt}'")
    return match.groupdict()


def read_coords(fid):
    """
    Read station coordinates from a comma- or whitespace-delimited text file.
    Blank lines and lines starting with '#' are ignored. If the first line
    contains a latitude column name, it is treated as a header and columns
    may be in any order (see `COORD_COLUMNS` for accepted names; unknown
    columns are ignored). Otherwise columns are assumed to be
    'station latitude longitude [elevation] [depth]'. e.g.,

        network,station,latitude,longitude,elevation
        XX,101,64.8591,-147.8492,140.5

    :type fid: str
    :param fid: path to the coordinate file
    :rtype: dict
    :return: {(network or None, station): (lat, lon, elevation, depth)};
        network is None if the file has no network column. Elevation [m]
        and depth [m] default to 0 if not given
    :raises ValueError: if required columns are missing or values are bad
    """
    aliases = {alias: name for name, names in COORD_COLUMNS.items()
               for alias in names}
    columns = ["station", "latitude", "longitude", "elevation", "depth"]

    with open(fid, "r") as f:
        lines = [line.strip() for line in f]
    lines = [line for line in lines if line and not line.startswith("#")]

    coords = {}
    for i, line in enumerate(lines):
        if "," in line:
            values = [val.strip() for val in line.split(",")]
        else:
            values = line.split()

        # Header row defines column order
        if i == 0 and any(aliases.get(val.lower()) == "latitude"
                          for val in values):
            columns = [aliases.get(val.lower()) for val in values]
            for name in ("station", "latitude", "longitude"):
                if name not in columns:
                    raise ValueError(f"'{fid}' header missing '{name}'")
            continue

        row = {name: val for name, val in zip(columns, values)
               if name and val}
        try:
            key = (row.get("network"), row["station"])
            coords[key] = (float(row["latitude"]), float(row["longitude"]),
                           float(row.get("elevation", 0)),
                           float(row.get("depth", 0)))
        except (KeyError, ValueError) as e:
            raise ValueError(f"'{fid}' could not parse line: '{line}'") from e

    return coords


def rename_stream(st, fid, fmt=None):
    """
    Set the trace stats (network, station, location, channel) of the first
    trace in a Stream from its file name. Does nothing if `fmt` is not given.

    :type st: obspy.core.stream.Stream
    :param st: stream to rename in place
    :type fid: str
    :param fid: file name the stream was read from
    :type fmt: str
    :param fmt: file name format string, see `parse_fname`
    """
    if not fmt:
        return
    parts = parse_fname(fid, fmt)
    for key in ("network", "station", "location", "channel"):
        if key in parts:
            setattr(st[0].stats, key, parts[key])


def build_inventory(get_template, fids, coords, fmt=None, sample_rate=None):
    """
    Build a comprehensive Inventory containing every unique channel found in
    a list of waveform files. Each channel receives a copy of the response
    of a template channel (e.g., from NRL) for its sample rate, coordinates
    from a user-defined coordinate file, an orientation from its channel
    code, and start and end dates spanning all data found for that channel.
    Channels recorded at more than one sample rate get one epoch per rate.

    :type get_template: callable
    :param get_template: function taking a sample rate [Hz] and returning
        the channel whose response, sensor and datalogger are copied
    :type fids: list of str
    :param fids: waveform file IDs, only the headers are read
    :type coords: dict
    :param coords: station coordinates, output of `read_coords`
    :type fmt: str
    :param fmt: optional file name format to set trace stats from file names,
        see `parse_fname`
    :type sample_rate: float
    :param sample_rate: expected sample rate [Hz] of the response; channels
        with a different sample rate are flagged with a warning
    :rtype: obspy.core.inventory.inventory.Inventory
    :return: inventory with one network per network code found in the data
    :raises ValueError: if any station is missing from `coords`
    """
    # Collect the time span of each unique channel and sample rate
    spans = {}
    for fid in fids:
        try:
            st = read(fid, headonly=True)
        except (TypeError, IsADirectoryError):
            print(f"{os.path.basename(fid)}... skipped, unknown file format")
            continue
        rename_stream(st, fid, fmt)
        for tr in st:
            key = (tr.id, tr.stats.sampling_rate)
            start, end = tr.stats.starttime, tr.stats.endtime
            if key in spans:
                start = min(start, spans[key][0])
                end = max(end, spans[key][1])
            spans[key] = (start, end)

    # Fail early if any station does not have coordinates
    missing = set()
    for tid, _ in spans:
        net, sta = tid.split(".")[:2]
        if (net, sta) not in coords and (None, sta) not in coords:
            missing.add(f"{net}.{sta}")
    if missing:
        raise ValueError(f"no coordinates for station(s): "
                         f"{', '.join(sorted(missing))}")

    # Build channels and sort them into {network: {station: [channels]}}
    channels = {}
    for (tid, rate), (start, end) in sorted(spans.items()):
        net, sta, loc, cha = tid.split(".")
        if sample_rate and rate != sample_rate:
            warnings.warn(f"{tid} sample rate {rate} Hz does not match "
                          f"response sample rate {sample_rate} Hz")
        try:
            template = get_template(rate)
        except ValueError as e:
            raise ValueError(f"{tid}: {e}") from e
        lat, lon, elev, depth = coords.get((net, sta),
                                           coords.get((None, sta)))
        azimuth, dip = ORIENTATIONS.get(cha[-1],
                                        (template.azimuth, template.dip))
        channel = Channel(code=cha, location_code=loc, latitude=lat,
                          longitude=lon, elevation=elev, depth=depth,
                          azimuth=azimuth, dip=dip, sample_rate=rate,
                          start_date=start, end_date=end,
                          sensor=deepcopy(template.sensor),
                          data_logger=deepcopy(template.data_logger),
                          response=deepcopy(template.response))
        channels.setdefault(net, {}).setdefault(sta, []).append(channel)

    networks = []
    for net, stations in channels.items():
        net_stations = []
        for sta, chans in stations.items():
            start = min(cha.start_date for cha in chans)
            end = max(cha.end_date for cha in chans)
            net_stations.append(
                Station(code=sta, latitude=chans[0].latitude,
                        longitude=chans[0].longitude,
                        elevation=chans[0].elevation, channels=chans,
                        start_date=start, end_date=end, creation_date=start)
                )
        networks.append(
            Network(code=net, stations=net_stations,
                    start_date=min(s.start_date for s in net_stations),
                    end_date=max(s.end_date for s in net_stations))
            )

    return Inventory(networks=networks, source="remresp.py",
                     created=UTCDateTime())


def build_fairfield_response(preamp_db=0, sample_rate=250, filter_phase="LP",
                             dc_filter="Off", output_units="mV", sensor_lf=5):
    """
    Build an NRL combine URL for a ZLand Gen2 sensor + datalogger cascade and
    return the resulting Inventory object with an attached response.

    Responses listed on NRL describe the generation 2 Zland Node datalogger. 
    Amplitudes are output from the datalogger in counts and scaled to 
    milliVolts by downloading software - responses are available for both 
    output unit types.

    This geophone measures ground velocity. NRL provides the onboard geophone 
    responses for the Zland Generation 2 Node instruments. They have published
    effective sensitivities and total damping values of 0.7, but resistances are 
    not specified. Zland Generation 1 Nodes used Geospace GS-30CT geophones 
    onboard.

    .. note ::

        There are two options for the Sensor response available on NRL. These
        are for a 10Hz corner and a 5Hz corner. The UAF sensors are 5Hz.
 
    :type preamp_db: int
    :param preamp_db: preamp gain [dB]; one of 0, 6, 12, 18, 24, 30, 36
    :type sample_rate: int
    :param sample_rate: final sample rate [Hz]; one of 250, 500, 1000, 2000
    :type filter_phase: str
    :param filter_phase: final FIR phase; 'LP' (linear) or 'MP' (minimum)
    :type dc_filter: str or int
    :param dc_filter: IIR low-cut setting; '1' (1 Hz) or 'Off'. Other
        corners (2-10 Hz) are not in the NRL; see NRL help for the pole edit
    :type output_units: str
    :param output_units: units of the data on disk; 'count' or 'mV'
    :type sensor_lf: int
    :param sensor_lf: geophone natural frequency [Hz]; 5 or 10
    :rtype: str
    :return: NRL combine URL returning cascaded StationXML
    """
    # Sensor sensivity options, will be hard-coded to 5Hz
    SENSOR_SENSITIVITY = {5: "76.7", 10: "78.7"}  # LF corner [Hz] -> V/(m/s)

    sensor = (f"sensor_MagseisFairfield_ZlandGen2Sensor_LF{sensor_lf}_"
              f"SG{SENSOR_SENSITIVITY[sensor_lf]}_STgroundVel")
    
    logger = (f"datalogger_MagseisFairfield_ZlandGen2_PD{preamp_db}_"
              f"FR{sample_rate}_FP{filter_phase}_DF{dc_filter}_"
              f"OU{output_units}")

    return (f"{BASE_NRL_URL}instconfig={sensor}:{logger}"
            f"&format=stationxml&nodata=404")


def build_smartsolo_response(preamp_db=0, sample_rate=250, filter_phase="LP",
                             dc_filter="Off", output_units="mV"):
    """
    Build an NRL combine URL for a SmartSolo IGU-BD3C-5

    The SmartSolo-IGU-BD3C-5 datalogger records three channels at a 
    preamplifier gain of 0 or 6 dB (gain factors 1 or 2) and a sample rate 
    of 50, 100, 125, 250, 500, 1000, 2000 or 4000 Hz. Amplitudes are output 
    from the datalogger in counts and scaled to milliVolts by downloading 
    software - responses are available for both output unit types. Its 
    onboard sensor is the three-component 5 second DT-Solo sensor. 

    DT-SOLO-BB: This intermediate period sensor measures ground velocity.
    It is the onboard sensors for the SmartSolo-IGU-BD3C-5.
 
    :type preamp_db: int
    :param preamp_db: preamp gain [dB]; one of 0, 6, 12, 18, 24, 30, 36
    :type sample_rate: int
    :param sample_rate: final sample rate [Hz]; one of 250, 500, 1000, 2000
    :type filter_phase: str
    :param filter_phase: final FIR phase; 'LP' (linear) or 'MP' (minimum)
    :type dc_filter: str or int
    :param dc_filter: IIR low-cut setting; '1' (1 Hz) or 'Off'. Other
        corners (2-10 Hz) are not in the NRL; see NRL help for the pole edit
    :type output_units: str
    :param output_units: units of the data on disk; 'count' or 'mV'
    :rtype: str
    :return: NRL combine URL returning cascaded StationXML
    """
    sensor = (f"sensor_DTCC_DT-SOLO-BB_LP5_SG209.4_STgroundVel")
    
    logger = (f"datalogger_DTCC_SmartSolo-IGU-BD3C-5_"
              f"PD{preamp_db}_FR{sample_rate}_FP{filter_phase}_"
              f"DF{dc_filter}_OU{output_units}")
    
    return (f"{BASE_NRL_URL}instconfig={sensor}:{logger}"
            f"&format=stationxml&nodata=404")


def get_nrl_inventory(args, sample_rate, cache=None):
    """
    Download the nominal response for the chosen instrument from NRL. If
    `args.sample_rate` is set it overrides `sample_rate`.

    :type args: argparse.Namespace
    :param args: parsed command line arguments, see `parse_args`
    :type sample_rate: float
    :param sample_rate: sample rate [Hz] of the data, e.g., from the header
    :type cache: dict
    :param cache: optional {sample_rate: Inventory} of previously downloaded
        responses, updated in place so each sample rate is fetched only once
    :rtype: obspy.core.inventory.inventory.Inventory
    :return: single-channel inventory with the attached response
    :raises ValueError: if NRL has no response for the sample rate
    """
    sample_rate = args.sample_rate or sample_rate
    if cache is not None and sample_rate in cache:
        return cache[sample_rate]

    # Header sample rates are floats, NRL expects one of the listed integers
    choices = PARAMETERS[args.choice]["sample_rate"]
    if sample_rate not in choices:
        raise ValueError(f"no NRL response for {sample_rate} Hz, must be one "
                         f"of {choices}")

    if args.choice == "fairfield":
        url = build_fairfield_response(preamp_db=args.pre_amp_gain,
                                       sample_rate=int(sample_rate),
                                       filter_phase=args.filter_phase,
                                       dc_filter=args.dc_filter,
                                       output_units=args.output_units,
                                       sensor_lf=5
                                       )
    elif args.choice == "smartsolo":
        url = build_smartsolo_response(preamp_db=args.pre_amp_gain,
                                       sample_rate=int(sample_rate),
                                       filter_phase=args.filter_phase,
                                       dc_filter=args.dc_filter,
                                       output_units=args.output_units
                                       )
    print(f"\n{url}\n")
    inv = read_inventory(url)
    print(f"{inv[0][0][0].response}\n")

    # Change dip to match orientation, see note above
    if args.choice == "fairfield":
        inv[0][0][0].dip = -90.0

    if cache is not None:
        cache[sample_rate] = inv

    return inv


def main():
    """
    Remove response from each file, using either a nominal NRL response or a
    user-provided StationXML file, or build a comprehensive StationXML file.
    """
    args = parse_args()

    assert(args.fids), f"{len(args.fids)} file IDs found"

    # Get response information from file, NRL responses are fetched per
    # sample rate as they are needed and cached by `get_nrl_inventory`
    from_nrl = args.choice != "stationxml"
    nrl_cache = {}
    if not from_nrl:
        inv = read_inventory(args.inv)
        print(f"\n{inv}\n")

    # Build and write one inventory for all files, no response removal
    if from_nrl and args.make_inv:
        coords = read_coords(args.coords)
        inv = build_inventory(
            lambda rate: get_nrl_inventory(args, rate, nrl_cache)[0][0][0],
            args.fids, coords, fmt=args.rename, sample_rate=args.sample_rate
            )
        inv.write(args.write_inv, format="STATIONXML")
        print(f"\n{inv}\n")
        print(f"Inventory written to '{args.write_inv}'")
        return

    if not os.path.exists(args.save):
        os.makedirs(args.save)

    # Begin response removal from each data stream
    print("\n")
    for i, fid in enumerate(args.fids):
        fid_out = os.path.basename(fid)
        print(fid_out, end="... ")

        # Check if this data has already been processed
        path_out = os.path.join(args.save, fid_out)
        if os.path.exists(path_out) and not args.overwrite:
            print("skipped, already processed")
            continue

        # Read data, ignore non waveform files
        try:
            st = read(fid)
        except (TypeError, IsADirectoryError):
            print("skipped, unknown file format")
            continue
    
        # Determine station naming from internal or from file
        rename_stream(st, fid, args.rename)

        # Rename the NRL inventory to match waveform so we can remove response
        if from_nrl:
            try:
                inv = get_nrl_inventory(args, st[0].stats.sampling_rate,
                                        nrl_cache)
            except ValueError as e:
                print(f"skipped, {e}")
                continue
            net, sta, loc, cha = st[0].id.split(".")
            inv[0].code = net
            inv[0][0].code = sta
            inv[0][0][0].code = cha
            inv[0][0][0].location_code = loc

        st.remove_response(inv, output=args.output, pre_filt=args.pre_filt,
                           water_level=args.water_level,
                           taper=bool(args.taper_fraction), 
                           taper_fraction=args.taper_fraction, 
                           zero_mean=not args.zero_mean_off)
        
        # Write out new file with response removed
        st.write(path_out, format="MSEED")
        print("complete")

        # Write Inventory and Response to file
        if i == 0 and from_nrl and args.write_inv:
            inv.write(args.write_inv, format="STATIONXML")
            print(f"Inventory written to '{args.write_inv}'")


if __name__ == "__main__":
    # See comment in top docstring about why this warning filter is here
    warnings.filterwarnings("ignore", message="The unit 'MV'")
    main()
