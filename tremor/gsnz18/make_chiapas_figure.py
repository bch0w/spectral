import os
import numpy as np
import matplotlib.pyplot as plt
from obspy import read, read_inventory, UTCDateTime, Stream


def read_write_around_chiapas(sta, plusminus=2):
    pathdir = ("/Users/chowbr/Documents/subduction/mseeds/BEACON/2017/XX/"
               "{sta}/HH{comp}.D/XX.{sta}.10.HH{comp}.D.2017.{jday}")
    chiapas = UTCDateTime("2017-09-08T00:00:00")
    inv = read_inventory("/Users/chowbr/Documents/subduction/mseeds/"
                         "BEACON/DATALESS/XX.RDF.DATALESS")
    for comp in ["N", "E"]:
        for jday in range(chiapas.julday-plusminus, chiapas.julday+plusminus):
            fid_in = pathdir.format(sta=sta, comp=comp, jday=jday)
            st = read(fid_in)
            st.attach_response(inv)
            st.remove_response(output="VEL", water_level=60, plot=False)
            fid_out = ("/Users/chowbr/Documents/subduction/spectral/tremor/"
                       "gsnz18/mseed_remove_response/{fid_in}")
            st.write(fid_out.format(
                fid_in=os.path.basename(fid_in)), format="mseed")


def read_around_chiapas(sta, plusminus=2):
    pathdir = ("/Users/chowbr/Documents/subduction/spectral/tremor/gsnz18/"
               "mseed_remove_response/XX.{sta}.10.HH{comp}.D.2017.{jday}")
    chiapas = UTCDateTime("2017-09-08T00:00:00")
    st = Stream()
    for comp in ["N", "E"]:
        for jday in range(chiapas.julday-plusminus, chiapas.julday+plusminus):
            fid_in = pathdir.format(sta=sta, comp=comp, jday=jday)
            st += read(fid_in)

    return st


def preprocess(st):
    st.merge()
    st.decimate(2)
    st.detrend("demean")
    st.detrend("linear")
    st.taper(max_percentage=0.05)

    return st


def combine_horizontals(st):
    st_north = st.select(component="N")
    st_east = st.select(component="E")
    st_out = st_north.copy()
    st_out.data = np.mean(st_north[0].data**2 + st_east[0].data**2)

    return st_out


def smoothed_envelope(st, window_length=10):
    samprate = st[0].stats.sampling_rate
    sample_window = int(window_length * samprate)
    rms_data = np.sqrt(st[0].data**2)
    envelope = []
    for i in range(0,len(rms_data)-sample_window,sample_window):
        window_average = sum(rms_data[i:i+sample_window])/sample_window
        envelope.append(window_average)
    envelope = np.array(envelope)

    return envelope


def amplitude_ratio(st, envelope=10, water_level=False):
    tremor_band = st.copy().filter('bandpass', freqmin=2, freqmax=5)
    quake_band = st.copy().filter('bandpass', freqmin=10, freqmax=20)
    surf_band = st.copy().filter('bandpass', freqmin=0.02, freqmax=0.1)
    if envelope:
        tremor_band = smoothed_envelope(tremor_band, window_length=envelope)
        quake_band = smoothed_envelope(quake_band, window_length=envelope)
        surf_band = smoothed_envelope(surf_band, window_length=envelope)
    else:
        tremor_band = tremor_band[0].data
        quake_band = quake_band[0].data
        surf_band = surf_band[0].data

    if water_level:
        def set_water_level(data):
            below_water_level = data < data.mean()
            data[below_water_level] = data.mean()

            return data
        quake_band = set_water_level(quake_band)
        surf_band = set_water_level(surf_band)

    ratio = tremor_band**2 / (quake_band * surf_band)
    two_sigma = np.std(ratio) * 2

    return ratio, two_sigma


def plot_waveforms():
    st = read_around_chiapas(sta='RD06')
    st = preprocess(st)
    st = combine_horizontals(st)
    ratio, two_sigma = amplitude_ratio(st, envelope=10)

    plt.plot(ratio)
    # plt.axhline(ratio, 0, 1)
    plt.show()

if __name__ == "__main__":
    plot_waveforms()








