"""
Remove response from Magseis-Fairfield ZLand 1C and 3C nodal data. 
Builds Node inventory with response information from Nominal Response Library.
Returns new MSEED files with the same name containing response removed data.

.. note:: Sign Convention

    Data converted with this script define the vertical axis as +Z up. Metadata
    should define `dip=-90` to maintain this orientation during response removal

    Fairfield Zland nodes assume a right handed coordinate system with +Z down
    (down positive) which is the standard in exploration seismology. 
    This is the opposite sense from earthquake/observational seismology, which 
    defines +Z up (up positive). ObsPy's read() function assumes the user is an 
    earthquake seismologist and flips (multiply by -1) the vertical axis so that 
    the output data are +Z up. To my understanding there is no accepted 
    convention for node users, and other groups may leave the vertical axis
    untouched, instead opting to set `dip=90` in the channel metadata to get
    data into the +Z up orientation.

.. rubric::
        
    python remresp.py <files> --output ./
"""
import os
import argparse
import numpy as np
from obspy import Inventory, UTCDateTime, read
from obspy.core.inventory.network import Network
from obspy.core.inventory.station import Station
from obspy.core.inventory.channel import Channel
from obspy.clients.nrl import NRL


ACCEPTABLE_SAMPLE_RATES = ["1000", "2000", "250", "500"]


def parse_args():
    """All modifications are accomplished with command line arguments"""
    parser = argparse.ArgumentParser()

    # Waveform Processing
    parser.add_argument("fids", nargs="+", help="required, file ID")
    parser.add_argument("-s", "--sample_rate", nargs="?", type=str, 
                        help="Sampling rate (Hz). If not given, sampling rate"
                             "will be taken from the first file, assuming that "
                             "all traces have the same sample rate",
                        choices=ACCEPTABLE_SAMPLE_RATES,
                        )
    parser.add_argument("-p", "--pre_amp_gain", nargs="?", type=str, 
                        default="18 dB (8)", 
                        help="Pre amplifier gain, typically '18 dB (8) but "
                             "talk to PI as it was set during job creation",
                        choices= ["0 dB (1)", "12 dB (4)", "18 dB (8)",
                                  "24 dB (16)", "30 dB (32)", "36 dB (64)", 
                                  "6 dB (2)"]
                        )
    parser.add_argument("-t", "--phase_type", nargs="?", type=str, 
                        default="Linear Phase", 
                        help="Final phase type. Typically 'Linear Phase'.", 
                        choices=["Linear Phase", "Minimum Phase"]
                        )
    parser.add_argument("-l", "--low_cut", nargs="?", type=str, 
                        default="Off", help="Low cut filter. Not usually used",
                        choices=["1 Hz", "Off"]
                        )
    parser.add_argument("-o", "--output", nargs="?", type=str, 
                        default="./response_removed",
                        help="where to save the newly created files")

    return parser.parse_args()


def return_response(pre_amp_gain="18 dB (8)", sample_rate="250", 
                    phase_type="Linear Phase", low_cut="Off"):
    """
    Build response from Nominal Response Library. Options for data logger which
    are Job dependent

    :type pre_amp_gain: str
    :param pre_amp_gain: '0 dB (1)', '12 dB (4)', '18 dB (8)', '24 dB (16)', 
                         '30 dB (32)', '36 dB (64)', '6 dB (2)'
    :type sample_rate: str
    :param sample_rate: '1000', '2000', '250', '500'
    :type phase_type: str
    :param phase_type: 'Linear Phase', 'Minimum Phase'
    :type low_cut: str
    :param low_cut: '1 Hz', 'Off'
    """
    nrl = NRL()
    sensor_keys = ["Magseis Fairfield", "Generation 2", "5 Hz"]
    datalogger_keys = ["Magseis Fairfield", "Zland 1C or 3C", pre_amp_gain, 
                       sample_rate, phase_type, low_cut]

    print(f"SENSOR KEYS: {sensor_keys}")
    print(f"DATALOGGER KEYS: {datalogger_keys}")

    response = nrl.get_response(sensor_keys=sensor_keys,
                                datalogger_keys=datalogger_keys)
    print(response)

    return response


def main():
    args = parse_args()

    # Few set up tasks
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    assert(args.fids), f"{len(fids)} file IDs found"

    if args.sample_rate is None:
        print(f"no sample rate given, retrieving from: {args.fids[0]}")
        st = read(args.fids[0])
        sample_rate = str(int(st[0].stats.sampling_rate))
        print(f"sample rate = {sample_rate}\n")
        assert(sample_rate in ACCEPTABLE_SAMPLE_RATES)
    else:
        sample_rate = args.sample_rate

    # Get response information
    response = return_response(pre_amp_gain=args.pre_amp_gain,
                               sample_rate=sample_rate,
                               phase_type=args.phase_type,
                               low_cut=args.low_cut) 

    # Dummy values to be used for Inventory, these are not actually important
    # for the response removal
    lat = 0.
    lon = 0.
    elevation = 0.
    depth = 0.
    start_date = UTCDateTime("1990-01-01")
    end_date = UTCDateTime()

    # Apply, remove and save
    for fid in args.fids:
        # Set up output file name and check existence
        fidout = os.path.basename(fid)
        pathout = os.path.join(args.output, fidout)
        if os.path.exists(pathout):
            print(f"{pathout} exists, skipping")
            continue
        else:
            print(fidout)

        st = read(fid)
        channels = []
        for tr in st:
            channel = Channel(code=tr.stats.channel, 
                              location_code=tr.stats.location, latitude=lat, 
                              longitude=lon, elevation=elevation, depth=depth,
                              response=response)
            channels.append(channel)
        station = Station(code=tr.stats.station, latitude=lat, longitude=lon,
                          elevation=elevation, channels=channels,
                          start_date=start_date, end_date=end_date
                          )
        network = Network(code=tr.stats.network, stations=[station]) 
        inv = Inventory(networks=[network])

        print("\tremoving response")
        st.remove_response(inv)

        print("\twriting file")
        st.write(pathout, format="MSEED")


if __name__ == "__main__":
    main()
