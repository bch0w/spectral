import os
import numpy as np
import matplotlib.pyplot as plt
from obspy import read, read_inventory, UTCDateTime, Stream
from matplotlib import gridspec

import sys
sys.path.append("../../modules")
from plotmod import pretty_grids, linespecs


def read_write_around_chiapas(sta, plusminus=2):
    """
    read in raw seismic data and response, remove response and save for later
    """
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
    """
    if read_write_around_chiapas() has been run, read in the preprocessed data
    """
    # pathdir = ("/Users/chowbr/Documents/subduction/spectral/tremor/gsnz18/"
    #            "mseed_remove_response/XX.{sta}.10.HH{comp}.D.2017.{jday}")
    pathdir = ("/seis/prj/fwi/bchow/spectral/tremor/gsnz18/"
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
    """
    amplitude envelope using RMS, smoothed out by averaging along window_length
    """
    samprate = st[0].stats.sampling_rate
    sample_window = int(window_length * samprate)
    rms_data = np.sqrt(st[0].data**2)
    envelope = []
    for i in range(0,len(rms_data)-sample_window,sample_window):
        window_average = sum(rms_data[i:i+sample_window])/sample_window
        envelope.append(window_average)
    envelope = np.array(envelope)

    return envelope


def amplitude_ratio(st, envelope=10, water_level=False, sigma_in=2):
    """
    tremor specific amplitude ratios based on frequency scanning method
    """
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
    sigma_out = np.std(ratio) * sigma_in

    return ratio, sigma_out


def values_over(values, over):
    """
    give a list of values over a certain threshold
    """
    under_over = values < over
    values[under_over] = np.nan

    return values


def plot_waveforms():
    st = read_around_chiapas(sta='RD06')
    st = preprocess(st)
    st = combine_horizontals(st)
    envelope_length = 60
    ratio, sigma = amplitude_ratio(st, sigma_in=2, envelope=envelope_length)
    t = np.linspace(0, (len(st[0])/st[0].stats.sampling_rate)/(3600),
                    len(ratio)
                    )
    f, ax = plt.subplots(1)
    linespecs()
    pretty_grids(ax)
    plt.plot(t, ratio, 'ko-', markersize=2)
    plt.plot(t, values_over(ratio, sigma), 'ro', markersize=2)
    plt.axhline(sigma, 0, 1)
    plt.show()


def mark_chiapas(origintime, time_break):
    chiapas = UTCDateTime("2017-09-08T04:49:00")
    mark = (chiapas - origintime)/time_break
    ax = plt.gca()
    ax.axvline(mark, 0, 1)


def gridspec_plot():
    """
    """
    st = read_around_chiapas(sta='RD06')
    st = preprocess(st)
    st = combine_horizontals(st)
    st_surfthewave = st.copy().filter('bandpass', freqmin=1/30, freqmax=1/6)
    envelope_length = 60
    ratio, sigma = amplitude_ratio(st, sigma_in=2, envelope=envelope_length)

    # plotting start
    f = plt.figure(figsize=(11.69, 8.27), dpi=75)
    linespecs()
    gs = gridspec.GridSpec(2, 1, height_ratios=[1,1], hspace=0)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharex=ax1)

    for ax in [ax1]:
        plt.setp(ax.get_xticklabels(), visible=False)

    time_break = 3600
    starttime = 0
    endtime = (len(st[0]) / st[0].stats.sampling_rate) / time_break

    pretty_grids(ax1)
    t = np.linspace(starttime, endtime, len(st[0].data))
    ax1.plot(t, st_surfthewave[0].data)

    pretty_grids(ax2)
    t = np.linspace(starttime, endtime, len(ratio))
    ax2.plot(t, ratio, 'ko-', markersize=2)
    ax2.plot(t, values_over(ratio, sigma), 'ro', markersize=2)
    ax2.axhline(sigma, 0, 1)
    mark_chiapas(st[0].stats.starttime, time_break)

    plt.show()



if __name__ == "__main__":
    gridspec_plot()








