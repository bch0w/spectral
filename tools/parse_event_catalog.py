"""
A script that is meant to parse through an ObsPy Catalog object, and pick out a
desired number of events that have a solid distribution of depths, magnitudes
and locations. Does this by taking statistical distributions of each variable
and choosing events based on these quantities
"""
import pandas as pd
from obspy import read_inventory
from obspy.geodetics.flinnengdahl import FlinnEngdahl


def get_event_id(event):
    """common method for deriving event id from event boject"""
    return event.resource_id.id.split('/')[1]


def parse_catalog_to_pandas(cat):
    """
    Turn the catalog object into a Pandas DataFrame
    :param cat:
    :return:
    """
    sources = pd.DataFrame()

    # Create a dataframe with source information, ignore duplicates
    for event in cat:
        event_id = get_event_id(event)
        latitude = event.preferred_origin().latitude
        longitude = event.preferred_origin().longitude
        nodal_planes = event.preferred_focal_mechanism().nodal_planes
        description = ""
        try:
            for dscrpt in event.event_descriptions:
                if dscrpt.type == "region name":
                    description = dscrpt.text.replace(", New Zealand", "")
                    break
            else:
                region = FlinnEngdahl().get_region(longitude, latitude)
                description = region.title().replace(", New Zealand", "")
        except IndexError:
            pass

        src = {
            "event_id": event_id,
            "time": str(event.preferred_origin().time),
            "magnitude": event.preferred_magnitude().mag,
            "depth_km": event.preferred_origin().depth * 1E-3,
            "latitude": latitude,
            "longitude": longitude,
            "location": description,
            "strike": nodal_planes.nodal_plane_1.strike,
            "dip": nodal_planes.nodal_plane_1.dip,
            "rake": nodal_planes.nodal_plane_1.rake,
            "strike_2": nodal_planes.nodal_plane_2.strike,
            "dip_2": nodal_planes.nodal_plane_2.dip,
            "rake_2": nodal_planes.nodal_plane_2.rake,
        }
        source = pd.DataFrame([list(src.values())],
                              columns=list(src.keys())
                              )
        source.set_index("event_id", inplace=True)

        sources = pd.concat([sources, source])

    return sources

