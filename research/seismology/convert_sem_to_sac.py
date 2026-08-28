"""
Convert SPECFEM ASCII synthetic waveforms (e.g., NN.SSS.CCC.semd) into SAC
files with appropriate SAC headers, using PySEP's `read_sem` utility.

If a CMTSOLUTION/SOURCE and STATIONS file are provided, SAC headers such as
event/station lat-lon, origin time and P/S arrival picks are also set.

Usage
-----
python convert_sem_to_sac.py --input_dir ./OUTPUT_FILES --output_dir ./sac \
    --source ./CMTSOLUTION --stations ./STATIONS
"""
import os
from glob import glob
from argparse import ArgumentParser

from pysep import read_sem, logger


def convert_sem_to_sac(input_dir, output_dir, glob_pattern="*.sem*",
                        source=None, stations=None, source_format="CMTSOLUTION",
                        origintime="1970-01-01T00:00:00"):
    """
    Read all SPECFEM ASCII synthetics matching `glob_pattern` in `input_dir`
    and write them out as SAC files of the same name in `output_dir`

    :type input_dir: str
    :param input_dir: directory containing SPECFEM ASCII synthetic files
    :type output_dir: str
    :param output_dir: directory to write out the resulting SAC files
    :type glob_pattern: str
    :param glob_pattern: wildcard pattern used to collect ASCII files,
        e.g., '*.sem*' catches '*.semd', '*.semv', '*.sem.ascii', etc.
    :type source: str
    :param source: path to SPECFEM source file (CMTSOLUTION, SOURCE, or
        FORCESOLUTION) used to set event-related SAC headers
    :type stations: str
    :param stations: path to SPECFEM STATIONS file used to set
        station-related SAC headers
    :type source_format: str
    :param source_format: format of `source`, e.g. 'CMTSOLUTION' or 'SOURCE'
    :type origintime: str
    :param origintime: fallback origin time used only if `source` is not
        given
    """
    os.makedirs(output_dir, exist_ok=True)

    fids = sorted(glob(os.path.join(input_dir, glob_pattern)))
    if not fids:
        logger.warning(f"no files matching '{glob_pattern}' found in "
                        f"{input_dir}")
        return

    for fid in fids:
        st = read_sem(fid, origintime=origintime, source=source,
                       stations=stations, source_format=source_format)

        fid_out = os.path.join(output_dir, f"{os.path.basename(fid)}.sac")
        st.write(fid_out, format="SAC")
        logger.info(f"{os.path.basename(fid)} -> {fid_out}")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input_dir", required=True,
                         help="directory containing SPECFEM ASCII synthetics")
    parser.add_argument("-o", "--output_dir", required=True,
                         help="output directory for converted SAC files")
    parser.add_argument("-g", "--glob_pattern", default="*.sem*",
                         help="wildcard pattern to match ASCII files "
                              "(default: '*.sem*')")
    parser.add_argument("-s", "--source", default=None,
                         help="path to CMTSOLUTION/SOURCE/FORCESOLUTION file "
                              "for event SAC headers")
    parser.add_argument("-f", "--source_format", default="CMTSOLUTION",
                         help="format of --source (default: 'CMTSOLUTION')")
    parser.add_argument("-t", "--stations", default=None,
                         help="path to STATIONS file for station SAC headers")
    parser.add_argument("--origintime", default="1970-01-01T00:00:00",
                         help="fallback origin time if --source not given")
    args = parser.parse_args()

    convert_sem_to_sac(
        input_dir=args.input_dir, output_dir=args.output_dir,
        glob_pattern=args.glob_pattern, source=args.source,
        stations=args.stations, source_format=args.source_format,
        origintime=args.origintime
    )
