#!/usr/bin/env python3
"""
Read in .fcnt data from the ZLand node and convert to MSEED data split on
full julian days with modified header values.

Data Format: https://docs.obspy.org/packages/obspy.io.rg16.html

.. note::

    This function assumes that nodes were oriented to geographic north, that
    the nodes are Fairfield ZLand nodes, and that each individual .fcnt file
    correspond to an individual instrument (two stations are NOT present in 
    the same file)

.. rubric::

    python fcnt2ms.py --files *.fcnt --network XX --output <path/to/output_dir>


Arguments:
    -f/--files <str>      # .fcnt files to read, wildcards okay
    -n/--network <str>    # two character network code for naming, defaults 'XX'
    -o/--output <str>     # output directory to save files, defaults to ./
    -c/--components <str> # optional list of components to select
    -b/--band_code <str>  # optional 'band' code if the guess function fails
    -i/--instrument_code <str>  # optional instrument code if not 'H'
"""
import argparse
import os
from obspy import read, UTCDateTime
from obspy.clients.nrl import NRL


def get_band_code(sampling_rate_hz):
    """
    Match sampling rate with SEED format channel code
    https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/

    .. note::

        Covers band codes L -> F and assumes we are using Fairfield nodes which
        have a corner frequency of 5 Hz

    :type dt: float
    :param dt: sampling rate of the data in seconds
    :rtype: str
    :return: band code as specified by Iris
    :raises KeyError: when dt is specified incorrectly
    """
    if sampling_rate_hz <= 1:  # 1Hz
        return "L"  # long period
    elif 1 < sampling_rate_hz < 10:  
        return "M"  # mid period
    elif 10 <= sampling_rate_hz < 80:  # [10, 80) Hz
        return "S"  # short period
    elif 80 <= sampling_rate_hz < 250:
        return "E"  # extremely short period
    elif 250 <= sampling_rate_hz < 1000:  # [250, 1000) Hz
        return "D"
    elif 1000 <= sampling_rate_hz < 5000:  # [1000, 5000) Hz
        return "G"
    else:
        raise KeyError("Channel code does not exist for this sampling rate, "
                       "please consult the IRIS SEED channel naming document "
                       "and input channel code through argument -b/--band_code"
                       )
    
def print_stream_info(st):
    """
    Print some information about a given Stream (read in from a single .fcnt 
    file) that may be useful for verifying the data is being read in correctly
    """
    # Determine all unique station names and components in this file 
    station_names = list(set(tr.stats.station for tr in st))

    # Get start and end times of the file to ensure we know what time 
    # range the data covered. Times will be in UTC
    tmin = min([tr.stats.starttime for tr in st])
    tmax = max([tr.stats.endtime for tr in st])
    julday_start = tmin.julday  # Julian day
    julday_end = tmax.julday  
    if tmin.year != tmax.year:
        raise Exception("Data spans multiple years, not supported, sorry!")

    # Print out some metadata for visual confirmation
    print(f"\t{len(st)} traces in Stream")
    print(f"\tstations: {station_names}")
    print(f"\tmin starttime: {tmin}")
    print(f"\tmax endttime:  {tmax}")
    print(f"\tjulian days range: {tmin.julday} -> {tmax.julday}")

def parse_st_for_jday_and_comp(st, jday, component):
    """
    Given a Stream object, return a new Stream object that contains only the
    traces that fall within the given julian day

    :type st: obspy.core.stream.Stream
    :param st: raw Stream object read in from .fcnt file
    :type jday: int
    :param jday: julian day to extract data for, e.g., 123
    :type component: str
    :param component: component to extract data for, e.g., 'Z'
    """
    st_out = st.select(component=component).copy()
    st_out.trim(starttime=UTCDateTime(f"{year}-{jday:0>3}T00:00:00Z"),
                endtime=UTCDateTime(f"{year}-{jday:0>3}T23:59:59.99999Z"),
                nearest_sample=False)
    
    st_out.merge()
    
    return st_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert .fcnt files to mseed")
    parser.add_argument("-f", "--files", type=str, nargs="+", 
                        help="List of .fcnt files to convert")
    parser.add_argument("-n", "--network", type=str, help="Network code to use",
                        default="XX")
    parser.add_argument("-c", "--components", type=str, default=None,
                         help="(OPTIONAL) only write out certain components, "
                              "e.g., 'ZN' will only write Z and N component. "
                              "If not provided, writes all available in data")
    parser.add_argument("-o", "--output", type=str, default=os.getcwd(),
                         help="output directory to save MSEED files. Defaults "
                              "to the current working directory")
    parser.add_argument("-O", "--overwrite", action="store_true", default=False,
                        help= "If data matching filename already exists, do "
                              "not overwrite it. If set to True, will "
                              "overwrite data")
    parser.add_argument("-b", "--band_code", type=str, default=None,
                        help="(OPTIONAL) override for the band code related to "
                             "sampling rate. this script includes a lookup "
                             "table but some codes require additional info "
                             "such as instrument corner")
    parser.add_argument("-i", "--instrument_code", type=str, default="H",
                        help="(OPTIONAL) override for the instrument code "
                             "related to the type of instrument. We assume "
                             "Fairfield nodes are 'H' for high gain "
                             "seismometer")

    # Get command line arguments ready
    args = parser.parse_args()
    fids = sorted(args.files)
    output = args.output
    overwrite = args.overwrite
    network = args.network
    components = args.components
    band_code = args.band_code
    instrument_code = args.instrument_code
    location = ""  # force to be empty but can be used if desired

    # Make sure the output directory exists
    os.makedirs(output, exist_ok=True)

    # Read in the .fcnt files
    print(f"READING {len(fids)} FILES FROM {fids[0]} -> {fids[-1]}")
    print(f"{'=' * 80}")
    for fid in fids:
        print(f"READING/WRITING FILE: {fid}")
        st = read(fid, format="rg16", contacts_north=True)
        print_stream_info(st)

        # Figure out what components we need to write files for
        if components:
            comp_list = list(components)
        else:
            comp_list = list(set(tr.stats.component for tr in st))

        # Figure out what julian days we need to write files for
        year = st[0].stats.starttime.year
        jday_start = st[0].stats.starttime.julday
        jday_end = st[-1].stats.endtime.julday
        
        # Small check to make sure we don't cross a year boundary
        year_stop = st[-1].stats.endtime.year
        if year != year_stop:
            raise Exception("Data spans multiple years, not supported, sorry!")

        # Sampling rate used to determine the band code, assuming same for all
        if band_code is not None:
            band_ = band_code
        else:
            band_ = get_band_code(st[0].stats.sampling_rate)
        
        station = st[0].stats.station

        for component in comp_list:
            channel = f"{band_}{instrument_code}{component}"  # e.g., GHZ
            for jday in range(jday_start, jday_end + 1):
                # Build the filename before doing any data manipulation so that
                # if we have already created the file we can skip right over
                filename = \
                   f"{network}.{station}.{location}.{channel}.{year}.{jday:0>3}"
                outfile = os.path.join(output, filename)
                if os.path.exists(outfile):
                    print(f"\tfile '{filename}' exists in 'output', skipping")
                    continue

                st_out = parse_st_for_jday_and_comp(st, jday=jday, 
                                                    component=component)
                for tr in st_out:
                    tr.stats.network = network
                    tr.stats.location = location  # drops default location code
                    tr.stats.channel = f"{band_}{instrument_code}{component}"

                print(f"\twriting file: {filename}")
                st_out.write(outfile, format="MSEED")
                
