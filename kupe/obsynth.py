import os
import sys
sys.path.append('../modules/')
from obspy import UTCDateTime, read, Stream

from getdata import pathnames, event_stream, get_moment_tensor
from synmod import stf_convolve, get_GCMT_solution
from procmod import preprocess
from plotmod import pretty_grids, align_yaxis

from os.path import join

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
print("===========================================")
print("==== Be aware: ignoring FutureWarnings ====")
print("===========================================")

# =================================== FUNC ====================================
def initial_data_gather(event_id,code):
    """gather event information, observation and synthetic traces,
    preprocess all traces accordingly and return one stream object with 6 traces
    """
    # station information
    network,station,location,channel = code.split('.')
    component = channel[-1]

    # event information
    MT = get_GCMT_solution(event_id)
    moment_tensor = MT.focal_mechanisms[0].moment_tensor
    half_duration = (moment_tensor.source_time_function['duration'])/2
    mt_geonet = get_moment_tensor(event_id)
    syn_start = UTCDateTime(mt_geonet['Date'])

    # grab synthetic data locally
    syntheticdata_path = join(pathnames()['syns'],event_id,'')
    syntheticdata = Stream()
    for c in ["N","E","Z"]:
        syntheticdata_filename = "{n}.{s}.BX{co}.semv.mseed".format(n=network,
                                                                s=station,
                                                                co=c)
        syntheticdata += read(join(syntheticdata_path,syntheticdata_filename))

    # grab observation data
    observationdata,inv,cat = event_stream(station=station,
                                            channel=channel,
                                            event_id=event_id,
                                            startpad=0,
                                            endpad=350)

    # preprocessing, instrument response, STF convolution
    observationdata_proc = preprocess(observationdata,inv)
    syntheticdata_preproc = preprocess(syntheticdata)
    syntheticdata_proc = stf_convolve(syntheticdata_preproc,half_duration)
    for tr in syntheticdata_proc:
        tr.stats.starttime = syn_start

    # combine, common sampling rate
    st_IDG = observationdata_proc + syntheticdata_proc
    st_IDG.resample(50)

    return st_IDG, MT


# =================================== MAIN ====================================
event_id = "2014p240655"
code = "NZ.PUZ..HHZ"
st,mt = initial_data_gather(event_id,code)

st.filter('bandpass',freqmin=1/30,freqmax=1/5)
st.plot()
