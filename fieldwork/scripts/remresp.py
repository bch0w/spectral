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
import argparse
import warnings
from obspy import read, read_inventory

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
    common.add_argument("--rename", action="store_true",
                        help="set NSLC stats from NN.SSSS.LL.CCC.YYYY.JJJ")
    common.add_argument("--write_inv", type=str, default="inventory.xml",
                        help="filename to write the response out as "
                             "a StationXML file for archival. This will write "
                             "the inventory for the first data file only. "
                             "Set as empty string to toggle off")

    # Instrument specific parameters
    parser = argparse.ArgumentParser(
        description="Remove nominal NRL response from node data")
    sub = parser.add_subparsers(dest="choice", required=True,
                                metavar="{fairfield,smartsolo}")

    helps = {"fairfield": "Magseis Fairfield ZLand 3C",
             "smartsolo": "SmartSolo IGU-BD3C-5"}
    for name, opts in PARAMETERS.items():
        p = sub.add_parser(name, parents=[common], help=helps[name])
        p.add_argument("-r", "--sample_rate", type=int, required=True,
                       choices=opts["sample_rate"],
                       help="final sample rate [Hz]")
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

    return parser.parse_args()


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


def main():
    args = parse_args()

    # Few set up tasks
    if not os.path.exists(args.save):
        os.makedirs(args.save)

    assert(args.fids), f"{len(args.fids)} file IDs found"

    # Get response information from NRL
    if args.choice == "fairfield":
        url = build_fairfield_response(preamp_db=args.pre_amp_gain, 
                                       sample_rate=args.sample_rate, 
                                       filter_phase=args.filter_phase,
                                       dc_filter=args.dc_filter, 
                                       output_units=args.output_units, 
                                       sensor_lf=5
                                       )
    elif args.choice == "smartsolo":
        url = build_smartsolo_response(preamp_db=args.pre_amp_gain, 
                                       sample_rate=args.sample_rate, 
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
        if args.rename:
            net, sta, loc, cha, *_ = fid.split(".")
            # Rename internal stats based on filename
            st[0].stats.network = net
            st[0].stats.station = sta
            st[0].stats.location = loc
            st[0].stats.channel = cha
        else:
            net, sta, loc, cha = st[0].id.split(".")

        # Rename the inventory to match waveform so we can remove response
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
        if i == 0 and args.write_inv:
            inv.write(args.write_inv, format="STATIONXML")
            print(f"Inventory written to '{args.write_inv}'")


if __name__ == "__main__":
    # See comment in top docstring about why this warning filter is here
    warnings.filterwarnings("ignore", message="The unit 'MV'")
    main()
