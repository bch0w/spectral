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
    -b/--band_code <str>  # optional 'band' code if the guess function fails
    -i/--instrument_code <str>  # optional instrument code if not 'H'
"""
import argparse
import os
from obspy import read, Stream, UTCDateTime


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert .fcnt files to mseed")
    parser.add_argument("-f", "--files", type=str, nargs="+", 
                        help="List of .fcnt files to convert")
    parser.add_argument("-n", "--network", type=str, help="Network code to use",
                        default="XX")
    parser.add_argument("-o", "--output", type=str, default=os.getcwd(),
                         help="output directory to save MSEED files. Defaults "
                              "to the current working directory")
    parser.add_argument("-O", "--overwrite", type=bool, default=False,
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
    band_code = args.band_code
    instrument_code = args.instrument_code

    # Make sure the output directory exists
    os.makedirs(output, exist_ok=True)

    # Read in the .fcnt files
    print(f"READING {len(fids)} FILES FROM {fids[0]} -> {fids[-1]}")
    print(f"{'=' * 80}")

    for fid in fids:
        print(f"READING/WRITING FILE: {fid}")
        st = read(fid, format="rg16", contacts_north=True)

        # Determine all unique station names and components in this file 
        station_names = list(set(tr.stats.station for tr in st))
        components = list(set(tr.stats.component for tr in st))

        # Print out some metadata for confirmation
        print(f"    {len(st)} traces in Stream")
        print(f"    stations: {station_names}")
        print(f"    components: {components}")
        print(f"    min starttime: {min([tr.stats.starttime for tr in st])}")
        print(f"    max endttime:  {max([tr.stats.endtime for tr in st])}")
        print(f"    begin parsing data...")

        # Sort by components, default are Z, N, E
        for component in components:
            st_comp = st.select(component=component)

            # Sort by time so we can just loop through and trim off days
            st_comp.sort(keys=["starttime"])

            # WARNING: I'm too lazy to figure out how to span a year, assuming 
            # we will not have deployments spanning New Years, which would 
            # increment the year value and break this code
            year = st_comp[0].stats.starttime.year
            year_end = st_comp[-1].stats.endtime.year
            if year != year_end:
                print("deployment spans a new year and this script is not "
                      "designed to handle, exiting")
                a = 1/0

            # Used to range over all available days
            jday_start = st_comp[0].stats.starttime.julday
            jday_end = st_comp[-1].stats.endtime.julday

            # Loop over days, extract data and write file for a given julian day
            jday = jday_start
            st_out = Stream()
            for i, tr in enumerate(st_comp):
                st_out.append(tr)
                if tr.stats.endtime.julday == jday:
                    continue
                # When we get to the Trace that starts on the next day, break
                else:
                    # Merge the files and trim on a full day
                    starttime = UTCDateTime(f"{year}-{jday:0>3}T00:00:00Z")
                    endtime = UTCDateTime(f"{year}-{jday:0>3}T23:59:59.99999Z")
                    st_out.merge()

                    # Trim away data that does not fall within this julian day
                    # do not use `nearest_sample` because that may pick up
                    # data from outside the specified bounds
                    st_out.trim(starttime, endtime, nearest_sample=False)

                    # Overwrite some metadata to be more SEED conforming
                    for tr_ in st_out:
                        # Change channel code to match standard SEED format
                        component = tr_.stats.component
                        if band_code is not None:
                            band_ = band_code
                        else:
                            band_ = get_band_code(tr_.stats.sampling_rate)

                        tr_.stats.network = network
                        tr_.stats.location = ""  # drop location code
                        tr_.stats.channel = \
                                f"{band_}{instrument_code}{component}"

                        # Generate file name based on metadata and timing
                        filename = (
                            f"{tr_.stats.network}.{tr_.stats.station}."
                            f"{tr_.stats.channel}.{tr_.stats.starttime.year}."
                            f"{tr_.stats.starttime.julday:0>3}"
                            )

                    # Optional overwriting command
                    if os.path.exists(filename) and not overwrite:
                        print(f"    file exists, won't overwrite: {filename}")
                        continue

                    # Write in MSEED Format
                    print(f"    writing mseed file: {filename}")
                    st_out.write(os.path.join(output, filename), format="MSEED")

                    # Reset stream and increment Julian day
                    st_out = Stream()

                    # If the current trace contains data from next day, we
                    # want to include it
                    if tr.stats.endtime.julday != jday:
                        st_out.append(tr)
                    
                    jday += 1


