"""
Read in .fcnt data from the ZLand node and convert to Mseed with proper metadata
https://docs.obspy.org/packages/obspy.io.rg16.html

.. note::

    This function assumes that nodes were oriented to geographic north 

.. rubric::

    python fcnt2ms.py --files *.fcnt --network XX 
"""
import argparse
from obspy import read, Stream, UTCDateTime


def channel_code(dt):
    """
    Match sampling rate with SEED format channel code
    https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/

    :type dt: float
    :param dt: sampling rate of the data in seconds
    :rtype: str
    :return: band code as specified by Iris
    :raises KeyError: when dt is specified incorrectly
    """
    if dt >= 1:
        return "L"  # long period
    elif 0.1 < dt < 1:
        return "M"  # mid period
    elif 0.0125 < dt <= 0.1:
        return "B"  # broad band
    elif 0.001 <= dt <= 0.0125:
        return "H"  # high broad band
    elif 0.004 <= dt < 0.001:
        return "C"
    elif 0.001 <= dt < 0.0002:
        return "F"
    else:
        raise KeyError("Channel code does not exist for this value of 'dt'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert .fcnt files to mseed")
    parser.add_argument("-f", "--files", type=str, nargs="+", 
                        help="List of .fcnt files to convert")
    parser.add_argument("-n", "--network", type=str, help="Network code to use",
                        default="XX")
    parser.add_argument("-s", "--station", type=str, default="Z",
                        help="Station prefix, e.g., 'Z' for Z01, Z02...")
    parser.add_argument("-c", "--components", type=str, default="ZNE")
    args = parser.parse_args()

    fids = args.files
    network = args.network
    station = args.station
    components = args.components

    # Read in the .fcnt files
    for fid in fids:
        print(f"reading {fid}")
        st = read(fid, format="rg16", contacts_north=True)

        # Sort by components, default are Z, N, E
        for component in components:
            st_comp = st.select(component=component)

            # Sort by time so we can just loop through and trim off days
            st_comp.sort(keys=["starttime"])

            # I'm too lazy to figure out how to span a year, assuming we will 
            # not have christmas + NY deployments
            year = st_comp[0].stats.starttime.year
            year_end = st_comp[-1].stats.endtime.year
            if year != year_end:
                print("deployment spans a new year and this script is not "
                      "designed to handle, exiting")
                a = 1/0

            # Used to range over all available days
            jday_start = st_comp[0].stats.starttime.julday
            jday_end = st_comp[-1].stats.endtime.julday

            # Loop over days, extract data and write, then delete data
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
                    st_out.trim(starttime, endtime)

                    # Overwrite some metadata to be more conforming with SEED
                    for tr_ in st_out:
                        tr_.stats.network = network
                        tr_.stats.station = f"{station}{tr_.stats.station[1:]}"
                        tr_.stats.location = ""

                        # Change channel code to match standard SEED format
                        # Assuming instr_ument code H for high gain seismometer
                        component = tr_.stats.component
                        tr_.stats.channel = \
                            channel_code(tr_.stats.delta) + "H" + component

                        # Generate file name based on metadata and timing
                        filename = (
                            f"{tr_.stats.network}.{tr_.stats.station}."
                            f"{tr_.stats.channel}.{tr_.stats.starttime.year}."
                            f"{tr_.stats.starttime.julday:0>3}"
                            )
                    st_out.write(filename, format="MSEED")
                    print(filename)

                    # Reset stream and increment julian day
                    st_out = Stream()

                    # If the current trace contains data from next day, we
                    # want to include it
                    if tr.stats.endtime.julday != jday:
                        st_out.append(tr)
                    
                    jday += 1

