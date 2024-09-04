"""
Generates rcv.txt and src.txt files for a nodal deployment, required by
the Fairfield HarvestManager software. Either take .kml files exported from
Google Earth or input lat/lon coordinates directly. 

.. note::

    To get .KML from Google Earth, create a directory, place points 
    (placemark) on the map geographically labelled with station names, and 
    then right click -> Save Place As... -> KML

.. rubric::

    python kml2srcrcv.py --file mar24_test.kml --line 1 --job_name mar24_test \
        --output rps
"""
import os
import argparse
import math
from pyproj import Proj


def utm_zone_from_lat_lon(lat, lon):
    """
    Calculate the UTM zone longitude value using quick maffs.
    Get the sign of the UTM zone based on the latitude value.

    .. note::

        Copied from Pyatoa.utils.srcrcv.utm_zone_from_lat_lon

    :type lat: float
    :param lat: latitude coordinate in degrees
    :type lon: float
    :param lon: longitude coordinate in degrees
    :rtype: int
    :return: UTM zone number
    """
    try:
        sign = lat / abs(lat)  # silly way to figure out if lat is +/-
    except ZeroDivisionError as e:
        raise Exception("latitude is 0, UTM zone is ambigious") from e
    return int(sign * math.ceil((lon + 180) / 6))


def latlon2utm(lat, lon):
    """
    Convert latitude and longitude coordinates to UTM projection using PyProj

    .. note::

        Modified from Pyatoa.utils.srcrcv.lonlat_utm

    :type lat: float or int
    :param lat: latitude value in WGS84 
    :type lon: float or int
    :param lon: longitude value in WGS84 
    :rtype: tuple (float, float, float)
    :return: (x in m UTM, y in m UTM, utm_zone)
    """
    # If converting latlon to utm and no utm zone given, calculate utm zone
    utm_zone = utm_zone_from_lat_lon(lat, lon)

    # Determine if the projection is north or south
    if utm_zone < 0:
        direction = "south"
    else:
        direction = "north"

    # Proj doesn't accept negative zones
    utm_zone_abs = abs(utm_zone)

    projstr = (f"+proj=utm +zone={utm_zone_abs} +{direction} +ellps=WGS84"
               f" +datum=WGS84 +units=m +no_defs")
    projection = Proj(projstr)

    x, y = projection(lon, lat, inverse=False)

    return x, y, utm_zone


def read_kml_file(fid):
    """
    Read a .kml file exported from Google Earth with point coordinates 
    specifying station locations. This function is super simple and just looks
    for line headers so it may not be very reliable.

    :type fid: str
    :param fid: name of the .kml file to read. See note for how to generate
    :rtype: dict of tuples
    :return: {<station_01>: (lat, lon)... }
    """
    stations = {}
    with open(fid, "r") as f:
        lines = f.readlines()

    record = False
    for line in lines:
        if "<Placemark" in line:
            record = True
        elif "</Placemark>" in line:
            record = False

        if record:
            if "<name>" in line:
                name = line.split("<name>")[1].split("</name>")[0]
            elif "<coordinates>" in line:
                lat = float(line.split("<coordinates>")[1].split(",")[1])
                lon = float(line.split("<coordinates>")[1].split(",")[0])
                stations[name] = (lat, lon)

    return stations


def write_src_rcv_txt(stations, job=None, rcv_fid_out="rcv.txt", 
                      src_fid_out="src.txt", line=1, path_out="./"):
    """
    Takes a Station file and writes out in format expected by Fairfield 
    HarvestManager software

    :type stations: dict of tuples
    :param stations: output of `read_kml_file`
    :type rcv_fid_out: str
    :param rcv_fid_out: name of the output receiver file, should be rcv.txt
    :type src_fid_out: str
    :param src_fid_out: name of output source file, should be src.txt
    :type line: int
    :param line: line number used for marking line and station
    """
    if job is not None:
        rcv_path_out = os.path.join(path_out, f"{job}_rcv.txt")
        src_path_out = os.path.join(path_out, f"{job}_src.txt")

    with open(rcv_path_out, "w") as f:
        f.write("Line Station UTMEast UTMNorth ElevMeters\n")  # header
        for i, (sta, coords) in enumerate(stations.items()):
            lat, lon = coords
            x_utm, y_utm, utm_zone = latlon2utm(lat, lon)
            print(f"{sta}: ({lat:.2f}, {lon:.2f}) -> "
                  f"({x_utm:.2f}, {y_utm:.2f}) UTM={utm_zone}")
            f.write(
                f"R{line} {line}{i:0>2} {x_utm:.2f} {y_utm:.2f} 0  # {sta}\n"
                )

    with open(src_path_out, "w") as f:
        f.write("Line Station UTMEast UTMNorth ElevMeters\n")  # header
        sta = list(stations.keys())[0]
        lat, lon = stations[sta]
        x_utm, y_utm, utm_zone = latlon2utm(lat, lon)
        f.write(f"S{line} {line}{i:0>2} {x_utm:.2f} {y_utm:.2f} 0\n")


def write_src_rcv_rps(stations, job=None, rcv_fid_out="rcv.rps", 
                      src_fid_out="src.rps", line=1, path_out="./"):
    """
    Write a tab-delimited SEG .rps file that can be used to import into
    the Fairfield software. See ZLAND quick start presentation for an example
    of what the file looks like. This method is preferred because it requires
    the least amount of mouse clicks to actually get the file to sync.

    :type stations: dict of tuples
    :param stations: output of `read_kml_file`
    :type rcv_fid_out: str
    :param rcv_fid_out: name of the output receiver file, should be rcv.txt
    :type src_fid_out: str
    :param src_fid_out: name of output source file, should be src.txt
    :type line: int
    :param line: line number used for marking line and station
    """
    line_format = "{prefix}{line:<16}{sta:<29}{x_utm:8.1f} {y_utm:8.1f} 0\n"

    if job is not None:
        rcv_path_out = os.path.join(path_out, f"{job}_{rcv_fid_out}")
        src_path_out = os.path.join(path_out, f"{job}_{src_fid_out}")
    else:
        rcv_path_out = os.path.join(path_out, rcv_fid_out)
        src_path_out = os.path.join(path_out, src_fid_out)
    
    with open(rcv_path_out, "w") as f:
        f.write("Line Station UTMEast UTMNorth ElevMeters\n")  # header
        for i, (sta, coords) in enumerate(stations.items()):
            lat, lon = coords
            x_utm, y_utm, utm_zone = latlon2utm(lat, lon)
            f.write(line_format.format(prefix="R", line=line, 
                                       sta=f"{line}{i:0>2}", x_utm=x_utm, 
                                       y_utm=y_utm)
                                       )

    with open(src_path_out, "w") as f:
        f.write("Line Station UTMEast UTMNorth ElevMeters\n")  # header
        sta = list(stations.keys())[0]
        lat, lon = stations[sta]
        x_utm, y_utm, utm_zone = latlon2utm(lat, lon)
        f.write(line_format.format(prefix="S", line=line, sta=f"{line}{i:0>2}", 
                                   x_utm=x_utm, y_utm=y_utm)
                                   )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert .kml to FF txt files")
    parser.add_argument("-f", "--file", type=str, nargs="?", 
                        help="List of .kml files to convert")
    parser.add_argument("-l", "--line", type=int, default=1,
                        help="Line number to assign to the given stations")
    parser.add_argument("-j", "--job_name", type=str, default=None,
                        help="Job name used to name the output files")
    parser.add_argument("-o", "--output", type=str, default="rps",
                        help="Output file format, choose one of 'rps' or 'txt'")
    args = parser.parse_args()

    stations = read_kml_file(args.file)
    if args.output == "txt":
        write_src_rcv_txt(stations, job=args.job_name, line=args.line)
    elif args.output == "rps":
        write_src_rcv_rps(stations, job=args.job_name, line=args.line)
    else:
        raise NotImplementedError(f"Output format {args.output} not recognized")
