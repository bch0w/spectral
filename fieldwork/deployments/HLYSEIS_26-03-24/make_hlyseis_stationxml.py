"""
Create a full StationXML with all necessary metadata for 2026-03 Healy Seismic
Deployment (HLYSEIS).

Station names and coordinates are read from the as-laid KMZ, and each station
is given a channel with the instrument response of the Magseis Fairfield
ZLand Gen2 node (5 Hz geophone, 18 dB preamp gain, 250 sps), taken from the
EarthScope Nominal Response Library (NRL).
"""
import os
import zipfile
import xml.etree.ElementTree as ET
from copy import deepcopy

from obspy import read_inventory, UTCDateTime
from obspy.core.inventory import (Inventory, Network, Station, Channel,
                                  Site, Equipment)


HERE = os.path.dirname(os.path.abspath(__file__))

# Deployment-specific parameters, edit as needed
KMZ_FILE = os.path.join(HERE, "2026-03-27_HLYSEIS_AS_LAID.kmz")
OUTPUT_FILE = os.path.join(HERE, "HLYSEIS_stationxml.xml")
NETWORK_CODE = "XX"  # temporary/unregistered network code
NETWORK_DESCRIPTION = "2026-03 Healy Seismic Deployment (HLYSEIS)"
SOURCE = "University of Alaska Fairbanks"
START_DATE = UTCDateTime("2026-03-25T00:00:00")
END_DATE = UTCDateTime("2026-04-25T00:00:00")  # pickup date
LOCATION_CODE = ""
SAMPLE_RATE = 250.  # must match the sample rate of the response

# Channel code: (azimuth, dip). Band code D = 250-1000 sps short period;
# instrument code P = geophone. Horizontals assume N oriented to true north.
CHANNELS = {"DPZ": (0., -90.), "DPN": (0., 0.), "DPE": (90., 0.)}

# Local copy of the NRL response; the NRL URL can also be used directly
RESPONSE_FILE = os.path.join(
    HERE, "..", "..", "metadata", "fairfield_PD18_FR250_FPLP_DFOff_OUmV.xml"
)
NRL_URL = (
    "https://service.earthscope.org/irisws/nrl/1/combine?instconfig="
    "sensor_MagseisFairfield_ZlandGen2Sensor_LF5_SG76.7_STgroundVel:"
    "datalogger_MagseisFairfield_ZlandGen2_PD18_FR250_FPLP_DFOff_OUmV"
    "&format=stationxml&nodata=404"
)
SENSOR_DESCRIPTION = "Magseis Fairfield ZLand Gen2 5 Hz geophone"
DATALOGGER_DESCRIPTION = "Magseis Fairfield ZLand Gen2, PD18, 250 sps"

# Station code prefix: site name, e.g., HC01 -> "Healy Creek 01"
SITE_NAMES = {"HC": "Healy Creek", "JB": "Jumbo Dome"}

KML_NS = {"kml": "http://www.opengis.net/kml/2.2"}


def read_kmz_stations(fid):
    """
    Read station names and coordinates from the Placemarks of a KMZ file.
    Coordinates are taken from each Placemark's <Point> element, which holds
    more precision than the rounded LAT/LON fields in the ExtendedData.

    :type fid: str
    :param fid: path to the KMZ (zipped KML) file
    :rtype: list of tuple
    :return: list of (name, latitude, longitude, elevation) for each station,
        with elevation in meters as given by the KML (often 0)
    """
    with zipfile.ZipFile(fid) as kmz:
        kml_name = [_ for _ in kmz.namelist() if _.endswith(".kml")][0]
        root = ET.fromstring(kmz.read(kml_name))

    stations = []
    for placemark in root.iterfind(".//kml:Placemark", KML_NS):
        name = placemark.findtext("kml:name", namespaces=KML_NS).strip()
        coords = placemark.findtext(".//kml:Point/kml:coordinates",
                                    namespaces=KML_NS)
        if coords is None:
            continue  # skip non-point features, e.g., lines or polygons
        # KML coordinate order is lon,lat[,alt]
        lon, lat, *alt = [float(_) for _ in coords.strip().split(",")]
        elv = alt[0] if alt else 0.
        stations.append((name, round(lat, 6), round(lon, 6), elv))

    names = [_[0] for _ in stations]
    assert len(names) == len(set(names)), f"duplicate station names in {fid}"

    return stations


def get_response(source=RESPONSE_FILE):
    """
    Read the instrument response from an NRL-generated StationXML file, which
    contains a single dummy channel (XX.YY.00.ZZZ) holding the response.

    :type source: str
    :param source: path or URL to the NRL StationXML file, defaults to the
        local copy in fieldwork/metadata/. Use `NRL_URL` to query the NRL
    :rtype: obspy.core.inventory.response.Response
    :return: the instrument response of the dummy channel
    """
    inv = read_inventory(source)
    channels = inv.get_contents()["channels"]
    assert len(channels) == 1, f"expected one channel in {source}"
    response = inv[0][0][0].response

    # Ensure the response matches the sample rate that we are setting
    response_sr = inv[0][0][0].sample_rate
    assert response_sr == SAMPLE_RATE, (
        f"response sample rate {response_sr} != {SAMPLE_RATE}"
    )

    return response


def get_site_name(code):
    """
    Expand a station code into a site name using the `SITE_NAMES` lookup,
    e.g., HC01 -> "Healy Creek 01"

    :type code: str
    :param code: station code made of a two-letter prefix and a number
    :rtype: str
    :return: expanded site name
    """
    prefix, number = code[:2], code[2:]
    assert prefix in SITE_NAMES, f"no site name for station prefix {prefix}"

    return f"{SITE_NAMES[prefix]} {number}"


def make_inventory(stations, response):
    """
    Build an ObsPy Inventory with one Station per entry in `stations`, and
    one Channel per entry in `CHANNELS`, each with a copy of `response`.

    :type stations: list of tuple
    :param stations: list of (name, latitude, longitude, elevation), as
        returned by `read_kmz_stations`
    :type response: obspy.core.inventory.response.Response
    :param response: instrument response to attach to each channel
    :rtype: obspy.core.inventory.inventory.Inventory
    :return: inventory containing all stations of the deployment
    """
    network = Network(code=NETWORK_CODE, description=NETWORK_DESCRIPTION,
                      start_date=START_DATE, end_date=END_DATE)

    for name, lat, lon, elv in stations:
        station = Station(code=name, latitude=lat, longitude=lon,
                          elevation=elv, start_date=START_DATE,
                          end_date=END_DATE, creation_date=START_DATE,
                          site=Site(name=get_site_name(name)))
        for cha_code, (azimuth, dip) in CHANNELS.items():
            channel = Channel(code=cha_code, location_code=LOCATION_CODE,
                              latitude=lat, longitude=lon, elevation=elv,
                              depth=0., azimuth=azimuth, dip=dip,
                              sample_rate=SAMPLE_RATE, start_date=START_DATE,
                              end_date=END_DATE)
            channel.sensor = Equipment(description=SENSOR_DESCRIPTION)
            channel.data_logger = Equipment(
                description=DATALOGGER_DESCRIPTION)
            channel.response = deepcopy(response)
            station.channels.append(channel)
        network.stations.append(station)

    return Inventory(networks=[network], source=SOURCE)


def main():
    """
    Read the as-laid KMZ, attach responses, and write StationXML to disk.
    """
    stations = read_kmz_stations(KMZ_FILE)
    print(f"read {len(stations)} stations from {os.path.basename(KMZ_FILE)}")

    response = get_response()
    inv = make_inventory(stations, response)
    print(inv)

    inv.write(OUTPUT_FILE, format="STATIONXML", validate=True)
    print(f"wrote {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
