"""
Remove response from Fairfield nodal data. Builds Node inventory with response 
information from Nominal Response Library
"""
import os
import argparse
import numpy as np
from obspy import Inventory, UTCDateTime, read
from obspy.core.inventory.network import Network
from obspy.core.inventory.station import Station
from obspy.core.inventory.channel import Channel
from obspy.clients.nrl import NRL


def parse_args():
    """All modifications are accomplished with command line arguments"""
    parser = argparse.ArgumentParser()

    # Waveform Processing
    parser.add_argument("fids", nargs="+", help="required, file ID")
    parser.add_argument("--pre_amp_gain", nargs="?", type=str, 
                        default="18 dB (8)", help="Pre amplifier gain")
    parser.add_argument("--sample_rate", nargs="?", type=str, default="250",
                        help="Sampling rate (Hz)")
    parser.add_argument("--phase_type", nargs="?", type=str, 
                        default="Linear Phase", help="Final phase type")
    parser.add_argument("--low_cut", nargs="?", type=str, 
                        default="Off", help="Low cut filter or not")
    parser.add_argument("--output", nargs="?", type=str, default="./remresp",
                        help="where to save new files")

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


if __name__ == "__main__":
    args = parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Get response information
    response = return_response(pre_amp_gain=args.pre_amp_gain,
                               sample_rate=args.sample_rate,
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

