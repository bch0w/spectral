"""
Read in .fcnt data from the ZLand node and convert to Mseed with proper metadata
https://docs.obspy.org/packages/obspy.io.rg16.html

.. note::

    This function assumes that nodes were oriented to geographic north 

.. rubric::

    python fcnt2ms.py --files *.fcnt --network XX 
"""
import argparse
import os
from glob import glob
from obspy import read, Stream, UTCDateTime


def channel_code(dt):
    """
    Match sampling rate with SEED format channel code

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


def read_data(fid, network=None, instrument_code="H"):
    """
    Read a single .fcnt file and return a formatted Stream object
    https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/

    :type fid: str
    :param fid: File path to .fcnt file
    :type network: str
    :param network: overwrite default network code. usually this is desired
        because the network codes will be the number of the line
    :type instrument_code: str
    :param instrument_code: instrument code as specified by IRIS, default is
        broadband high gain seismometer 'H'
    :rtype: obspy.core.stream.Stream
    :return: Stream object with metadata and data filled in
    """
    st = read(fid, format="rg16", contacts_north=True)
    st.merge()  # merge 30-second chunks into continous waveforms

    for tr in st:
        # overwrite network code
        if network:
            tr.stats.network = network
        # change channel code to match standard SEED formatting
        dt = tr.stats.delta
        component = tr.stats.component
        # e.g., HHZ, HHN, HHE
        tr.stats.channel = channel_code(dt) + instrument_code + component

        # strip location code, not required
        tr.stats.location = ""

    return st
        

def split_st_on_jdays(st):
    """
    Split a continuous Stream into day long sections for writing to mseed files
    
    .. warning::

        No safeguard for deployments that cross over New Years. I'm assuming
        were not out in the field at that time anyway
    """
    year = st[0].stats.starttime.year
    jday_start = st[0].stats.starttime.julday
    jday_end = st[0].stats.endtime.julday

    st_out = Stream()
    for jday in range(jday_start, jday_end + 1):
        starttime = UTCDateTime(f"{year}-{jday:0>3}T00:00:00Z")
        endtime = UTCDateTime(f"{year}-{jday+1:0>3}T00:00:00Z")
        st_out.extend(st.slice(starttime, endtime))

    return st_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert .fcnt files to mseed")
    parser.add_argument("-f", "--files", type=str, nargs="+", 
                        help="List of .fcnt files to convert")
    parser.add_argument("-n", "--network", type=str, help="Network code to use")
    args = parser.parse_args()

    fids = args.files
    network = args.network

    # Read in the .fcnt files
    for fid in fids:
        print(fid)
        st = read_data(fid, network=network)
        st = split_st_on_jdays(st)

        # Write out the Stream to mseed files
        for tr in st:
            # Generate file name based on metatdata and timing
            filename = (f"{tr.stats.network}.{tr.stats.station}."
                        f"{tr.stats.channel}.{tr.stats.starttime.year}."
                        f"{tr.stats.starttime.julday:0>3}")
            tr.write(filename, format="MSEED")
