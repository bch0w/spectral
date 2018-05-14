"""to build up a catalog of earthquakes useable for seismic tomography

potential workflow outline
-collect all stations on the north island, as well as temporary deployments,
offshore deployments, rdf etc.
-for each station, collect stationxml information into a single xml file, and
numpy arrays containing lat/lon of each station
-collect all m4 or m5 events in the area; function ability to manually look at
waveforms and assess quality
-collect all <m3 events for a radius surrounding each station for more events,
using lat/lon as inputs to event getter
-collect all quakeml information into single quakeml file, somehow organize
which events belong to which stations (excel file?)
-run detection algorithm on all station data?
-run flexwin on all picks for optimum window choice?
"""
import sys
sys.path.append('../modules')
import traceback
import collections
import pandas as pd
import numpy as np
from obspy import UTCDateTime, read_events
from obspy.clients.fdsn import Client

import getdata
import synmod
from getdata import pathnames

def collect_events():
    c = Client("GEONET")
    new_zealand = [-50,-35,165,180]
    kaikoura_to_east_cape = [-43,-37,172,180]
    blenheim_to_east_cape = [-41.5,-37,172,180]
    lat_lon = kaikoura_to_east_cape
    cat = c.get_events(starttime="2005-01-01T00:00:00",
                        endtime=UTCDateTime(),
                        minmagnitude=4.5,
                        maxmagnitude=6,
                        maxdepth=80,
                        minlatitude=lat_lon[0],
                        maxlatitude=lat_lon[1],
                        minlongitude=lat_lon[2],
                        maxlongitude=lat_lon[3],
                        orderby="magnitude")
    return cat

def info_from_GCMT(event_id):
    """modified from synmod.tshift_halfdur()
    get centroid time, half duration and GCMT moment tensor solution
    """
    GCMT = getdata.get_GCMT_solution(event_id)
    # time
    centroid_time = [i.time for i in GCMT.origins
                                            if i.origin_type == "centroid"][0]
    # magnitude
    mwc = [i.mag for i in GCMT.magnitudes
                                        if i.magnitude_type == "Mwc"][0]

    # half duration
    moment_tensor = GCMT.focal_mechanisms[0].moment_tensor
    half_duration = (moment_tensor.source_time_function['duration'])/2

    # GCMT moment tensor, convert to dyne*cm
    GCMT_RTP = synmod.mt_from_event(GCMT)
    for entry in GCMT_RTP:
        GCMT_RTP[entry] = GCMT_RTP[entry]*1E7

    return centroid_time, half_duration, GCMT_RTP, mwc

def info_from_RISTAU(event_id):
    """taken from generate_CMTSOLUTION
    get information from Ristau's MT solution
    """
    MT = getdata.get_moment_tensor(event_id)
    mt = [MT['Mxx'],MT['Myy'],MT['Mzz'],MT['Mxy'],MT['Mxz'],MT['Myz']]
    mt = [_*(1E20) for _ in mt]
    mt = synmod.mt_transform(mt,method='xyz2rtp')
    mrr,mtt,mpp,mrt,mrp,mtp = mt
    GEONET_MT_DICT = collections.OrderedDict({"m_rr":mrr,
                                              "m_tt":mtt,
                                              "m_pp":mpp,
                                              "m_rt":mrt,
                                              "m_rp":mrp,
                                              "m_tp":mtp})

    mw = MT['Mw']

    return GEONET_MT_DICT, mw

def parse_catalog(cat):
    """takes an obspy catalog object and parses out all event information into
    relevant lists, included information useful for tomographic simulations:
    + event_id: taken from GEONET event identification
    + date: in UTCDateTime format
    + depth: in units of km
    + lat/lon: included with uncertainties just incase
    - RTP MT: Ristau XYZ MT converted to GCMT coordinate system [dyne*cm]
        NOTE: will have the labelling m_??
    - GCMT: RTP MT from GCMT, for comparisons with converted Ristau solution
        NOTE: will have the labelling M??
    - centroid time: from GCMT
    - half duration: from GCMT
    """
    event_ids,dates,lats,lat_uncs,lons,lon_uncs,depths = [],[],[],[],[],[],[]
    tshifts,centroid_times,half_durations,GCMT_RTPs,GEONET_RTPs = [],[],[],[],[]
    mws,mwcs,delta_m = [],[],[]
    errors,exceptions = [],[]

    # NaN dictionaries
    GEONET_NaN = collections.OrderedDict({"m_rr":np.nan,
                                          "m_tt":np.nan,
                                          "m_pp":np.nan,
                                          "m_rt":np.nan,
                                          "m_rp":np.nan,
                                          "m_tp":np.nan
                                          })
    GCMT_NaN = collections.OrderedDict({"Mrr":np.nan,
                                        "Mtt":np.nan,
                                        "Mpp":np.nan,
                                        "Mrt":np.nan,
                                        "Mrp":np.nan,
                                        "Mtp":np.nan
                                        })


    for event in cat:
        # information from GEONET FDSN catalog
        event_id = str(event.resource_id).split('/')[1]
        print(event_id)
        try:
            origins = event.origins[0]
            date = origins.time
            depth = origins.depth * 1E-3
            lat = origins.latitude
            lon = origins.longitude
            lat_uncertainty = origins.latitude_errors['uncertainty']
            lon_uncertainty = origins.longitude_errors['uncertainty']
        except Exception as e:
            errors.append(event_id)
            exceptions.append(e)
            print('\t> error EVENT')
            continue


        # information from GEONET (Ristau) solutions
        try:
            strike1 = False
            GEONET_RTP, mw = info_from_RISTAU(event_id)
        except Exception as e:
            errors.append(event_id)
            exceptions.append(e)
            print('\t> error GEONET')
            GEONET_RTP = GEONET_NaN
            mw = np.nan
            strike1 = True
            pass

        # information from GCMT solution
        try:
            centroid_time, half_duration, GCMT_RTP, mwc = info_from_GCMT(
                                                                    event_id)
        except Exception as e:
            if strike1:
                continue
            errors.append(event_id)
            exceptions.append(e)
            print('\t> error GCMT')
            centroid_time = np.nan
            half_duration = np.nan
            GCMT_RTP = GCMT_NaN
            mwc = np.nan
            pass

        # try to create some tshift entry
        try:
            tshift = abs(date-centroid_time)
        except Exception as e:
            tshift = np.nan
            pass

        # append to lists
        event_ids.append(event_id)
        dates.append(date)
        lats.append(lat)
        lons.append(lon)
        lat_uncs.append(lat_uncertainty)
        lon_uncs.append(lon_uncertainty)
        depths.append(depth)

        mws.append(mw)
        mwcs.append(mwc)
        delta_m.append(abs(mw-mwc))
        tshifts.append(tshift)

        GEONET_RTPs.append(GEONET_RTP)
        GCMT_RTPs.append(GCMT_RTP)
        centroid_times.append(centroid_time)
        half_durations.append(half_duration)

    # errors
    errorCat = pd.DataFrame({"event_id" : errors,
                             "exception" : exceptions
                             })

    # initialize tomCat
    baseCat = pd.DataFrame(collections.OrderedDict({
                                        "event_id" : event_ids,
                                        "datetime" : dates,
                                        "cetroid_time" : centroid_times,
                                        "mw" : mws,
                                        "mwc" : mwcs,
                                        "delta_m" : delta_m,
                                        "tshift" : tshifts,
                                        "half_dur" : half_durations,
                                        "latitude" : lats,
                                        "longitude" : lons,
                                        "d_lat" : lat_uncs,
                                        "d_lon" : lon_uncs,
                                        "depth" : depths,
                                        })
                                        )

    # moment tensor dataframes
    geonetRTPs = pd.DataFrame(GEONET_RTPs)
    gcmtRTPs = pd.DataFrame(GCMT_RTPs)
    RTPs = geonetRTPs.join(gcmtRTPs)
    tomCat = baseCat.join(RTPs)

    tomCat = tomCat.sort_values(by="datetime")
    errorCat = errorCat.sort_values(by="event_id")

    return tomCat,errorCat


if __name__ == "__main__":
    catpath = pathnames()['data'] + 'QUAKEML/kaikoura_to_east_cape_update.xml'
    # catpath = pathnames()['data'] + 'QUAKEML/catbuild_testcat.xml'
    cat = read_events(catpath)
    tomCat,errorCat = parse_catalog(cat)

    # save as pickle and csv files
    outpath = pathnames()['kupe'] + 'tomCat/{}{}'
    tomCat.to_pickle(outpath.format("tomCat",""))
    errorCat.to_pickle(outpath.format("errorCat",""))
    tomCat.to_csv(outpath.format("tomCat",".csv"),index=False)
    errorCat.to_csv(outpath.format("errorCat",".csv"),index=False)
