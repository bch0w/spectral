import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from obspy import read, read_inventory, UTCDateTime, Stream
from matplotlib import gridspec
import matplotlib as mpl

import sys
sys.path.append("../../modules")
from plotmod import pretty_grids, linespecs
sys.path.append("..")
from utils import create_min_max
from obspy.signal.cross_correlation import correlate

mpl.rcParams['font.size'] = 15
mpl.rcParams['lines.linewidth'] = 1.
mpl.rcParams['lines.markersize'] = 1.75
mpl.rcParams['axes.linewidth'] = 2.5

def detect_earthquakes(st,corr_criteria=0.7):
    """try to remove earthquake from waveforms by taking correlations with
    an exponential function. If correlation criteria met, earthquake 'detected'
    :type de_array: numpy array
    :param de_array: datastream representing waveform envelope
    :type sampling_rate: float
    :param sampling_rate: sampling rate
    :type corr_criteria: float
    :param corr_criteria: threshold for detecting earthquakes, defaults to 0.7
    :rtype quakearray: np.array
    :return quakearray: de_array containing -1's for detected earthquakes
    """
    # set exponential template
    sampling_rate = int(st[0].stats.sampling_rate)
    de_array = st[0].data

    sampling_rate_min = sampling_rate * 60
    sampling_rate_half_min = int(sampling_rate_min * (1/2))
    sampling_rate_one_one_half_min = int(sampling_rate_min * (3/2))
    x= np.linspace(0.002,6,sampling_rate_min)
    exp_internal = -(x/2)*2
    exp_template = np.exp(exp_internal)

    # fill-value arrays if earthquake detected
    nan_fill = np.nan*(np.ones(sampling_rate_min))
    nan_fill_ext = np.nan*(np.ones(sampling_rate_one_one_half_min))

    quakecount = 0
    quakearray = np.array([])
    for S0 in range(0,len(de_array),sampling_rate_min):
        S1 = S0 + sampling_rate_min
        tremor_snippet = de_array[S0:S1]
        exp_correlation = correlate(a=exp_template,
                                    b=tremor_snippet,
                                    shift=len(x))

        if exp_correlation.max() > corr_criteria:
            quakecount +=1
            if S0 == 0:
                quakearray = np.append(quakearray,nan_fill)
            else:
                quakearray_new = quakearray[:S0-sampling_rate_half_min]
                quakearray = np.append(quakearray_new,nan_fill_ext)
        else:
            quakearray = np.append(quakearray,tremor_snippet)

    return quakearray

def read_around_chiapas(sta, plusminus=2):
    """
    if read_write_around_chiapas() has been run, read in the preprocessed data
    """
    pathdir = ("/Users/chowbr/Documents/subduction/mseeds/CHIAPAS/"
               "??.{sta}.10.??{comp}.D.2017.{jday}")
    # pathdir = ("/seis/prj/fwi/bchow/spectral/tremor/gsnz18/"
    #            "mseed_remove_response/XX.{sta}.10.HH{comp}.D.2017.{jday}")
    chiapas = UTCDateTime("2017-09-08T00:00:00")
    st = Stream()
    for comp in ["N", "E"]:
        for jday in range(chiapas.julday-plusminus, chiapas.julday+plusminus):
            fid_in = glob.glob(pathdir.format(sta=sta, comp=comp, jday=jday))
            if fid_in:
                st += read(fid_in[0])

    return st


def preprocess(st):
    st.merge(fill_value=0)
    st.decimate(2)
    st.detrend("demean")
    st.detrend("linear")
    st.taper(max_percentage=0.05)

    return st


def combine_horizontals(st):
    if len(st[0]) != len(st[1]):
        starttime = UTCDateTime("{yr}-{jday}T00:00:00".format(
            yr=st[0].stats.starttime.year, jday=st[0].stats.starttime.julday+1)
        )
        endtime = UTCDateTime("{yr}-{jday}T00:00:00".format(
            yr=st[0].stats.endtime.year, jday=st[0].stats.endtime.julday)
        )
        st.trim(starttime=starttime, endtime=endtime)
    st_north = st.copy().select(component="N")
    st_east = st.copy().select(component="E")
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
    for i in range(0, len(rms_data)-sample_window, sample_window):
        # window_average = sum(rms_data[i:i+sample_window])/sample_window
        window_median = np.median(rms_data[i:i+sample_window])
        envelope.append(window_median)
    envelope = np.array(envelope)

    return envelope


def amplitude_ratio(st, envelope=10, water_level=False, sigma_in=2):
    """
    tremor specific amplitude ratios based on frequency scanning method
    """
    tremor_band = st.copy().filter('bandpass', freqmin=2, freqmax=8)
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


def filter_save():
    for sta in ["ANWZ", "PRHZ", "RD01", "RD02", "RD03", "RD05", "RD16", "BFZ"]:
        try:
            st = read_around_chiapas(sta=sta, plusminus=3)
            st = preprocess(st)
            # 2-5 Hz data
            st_filtered = st.copy().filter('bandpass', freqmin=2,
                                           freqmax=8
                                           )
            for comp in ["N", "E"]:
                st_out = st_filtered.select(component=comp)
                fid_out = ("/Users/chowbr/Documents/subduction/mseeds/CHIAPAS/"
                           "filtered_2-8hz/{}".format(st_out[0].get_id())
                           )
                st_out.write(fid_out, format='MSEED')
        except Exception as e:
            print(e)
            continue

        #     st = combine_horizontals(st)
        #     x, y = create_min_max(st_filtered[0])
        #     y /= y.max()
        #     y += step
        #     plt.plot(x,y)
        #     step += 1


def time_to_decimal(datetime):
    return datetime.julday + (datetime.hour + datetime.minute/60)/24


def plot_ratios():
    f, ax = plt.subplots(1)
    sta = "RD02"
    # pathdir = ("/Users/chowbr/Documents/subduction/mseeds/CHIAPAS/"
    #            "??.{sta}.10.??{comp}.D.2017.{jday}")
    pathdir = ("/seis/prj/fwi/bchow/mseeds/CHIAPAS/XX.{sta}.10.HH{comp}.D.2017.{jday}")
    chiapas = UTCDateTime("2017-09-08T00:00:00")
    quiet = UTCDateTime("2017-12-09T00:00:00")
    for day, color in zip([chiapas, quiet],['r','k']):
        st = Stream()
        for comp in ["N", "E"]:
            # import ipdb;ipdb.set_trace()
            fid_in = glob.glob(pathdir.format(sta=sta, comp=comp,
                                              jday=day.julday))
            if fid_in:
                st += read(fid_in[0])
        st = preprocess(st)
        st = combine_horizontals(st)
        envelope_length = 5 * 60
        ratio, sigma = amplitude_ratio(st, sigma_in=2,
                                       envelope=envelope_length,
                                       water_level=True
                                       )
        # x = np.linspace(time_to_decimal(st[0].stats.starttime),
        #                 time_to_decimal(st[0].stats.endtime),
        #                 len(ratio)
        #                 )
        x = np.linspace(0,24,len(ratio))
        plt.plot(x, ratio, '{}o-'.format(color), markersize=2,
                 label="{}.{}".format(day.year, day.julday), markerfacecolor=color,
                 markeredgewidth=1.2, linewidth=1.5)
        plt.plot(x, values_over(ratio, sigma), 'ko', markerfacecolor='w',
                 markersize=5, markeredgewidth=1.2)

            # except Exception as e:
            #     print(e)
            #     continue
    # plt.axvline(chiapas, 0, 1, c='r', zorder=2)
    pretty_grids(ax)
    plt.legend()
    plt.xlabel('Hours')
    plt.xlim([0,24])
    plt.ylabel('Amplitude Ratio R(t)')
    # plt.setp(ax.get_yticklabels(), visible=False)

    plt.show()


def gridspec_plot():
    """
    """
    st = read_around_chiapas(sta='PRHZ', plusminus=3)
    # import ipdb;ipdb.set_trace()

    st = preprocess(st)
    st = combine_horizontals(st)
    st_surfthewave = st.copy().filter('bandpass', freqmin=2, freqmax=5)
    envelope_length = 60*60
    ratio, sigma = amplitude_ratio(st, sigma_in=2, envelope=envelope_length)
    sigma = 0.0675 + 0.077

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
    plot_ratios()








