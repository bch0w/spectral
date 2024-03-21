"""
Caclulate the minimum resolvable period of a numerical mesh by comparing
synthetic waveforms generated using a coarse (ngll5) and fine (ngll7) mesh.
"""
import os
import traceback
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from glob import glob
from obspy import Stream, UTCDateTime
from obspy.signal.cross_correlation import correlate, xcorr_max

from pyatoa import read_sem
from pyatoa.visuals.wave_maker import pretty_grids


dummy_time = UTCDateTime("2000-01-01T00:00:00")


class Minresper:
    """
    The Minimum Resolvable Period class.
    """
    def __init__(self, path_ngll5, path_ngll7):
        """
        Set internally used parameters such as path and station names.
        """
        self.path_ngll5 = path_ngll5
        self.path_ngll7 = path_ngll7
        self.stations = self.get_station_wildcards()
        self.stats = {"code": [], "comp": [], "period": [], "value": []}
        self.dfstats = None
        self._current = None

        # For plotting
        self.f = None
        self.ax = None

    def get_station_wildcards(self):
        """
        Get a list of stations with wildcards in the channels so that the
        different sampling rates between NGLL5 and NGLL7 can be addressed using
        glob.

        :rtype: np.array
        :return: wildcard station names
        """
        codes = [os.path.splitext(os.path.basename(_))[0] for _ in
                 glob(os.path.join(self.path_ngll5, "*sem?"))
                 ]

        # Replace the channel codes with wildcard strings
        for i, code in enumerate(codes):
            cha = code.split(".")[-1]
            codes[i] = code.replace(cha, "?X?")

        return np.unique(codes)

    def initiate_plot(self, dpi=100):
        """
        Simple setup for plotting routines.
        """
        f, ax = plt.subplots(figsize=(500/dpi, 800/dpi), dpi=dpi)
        self.f = f
        pretty_grids(ax)
        self.ax = ax

    def finalize_plot(self, show=True, save=None, **kwargs):
        """
        Final touches to the plot before showing or saving.
        """
        plt.sca(self.ax)
        # plt.title(f"{self._current} Minimum Resolvable Period")
        plt.title(f"{self._current} Model Interpolation Differences")
        plt.grid()
        plt.legend(loc="upper right")
        plt.xlabel("Time [s]")
        # plt.ylabel("Filter bands [s], MRP Value")
        plt.ylabel("Filter bands [s], Max CC Value")
        plt.gca().axes.get_yaxis().set_ticks([])

        if save:
            plt.savefig(os.path.join(save, f"{self._current}.png"))
        if show:
            plt.show()
        else:
            plt.close()

    def plot_traces(self, tr_a, tr_b, step=0, c_ngll5="k", c_ngll7="r", tag=""):
        """
        Convenience function to plot the two NORMALIZED traces with a step added
        to the data so that multiple traces can be plotted together.
        """
        tr_a_plot = tr_a.copy()
        tr_b_plot = tr_b.copy()

        # Only need to label once for legend
        if step == 0:
            label_a = "Coarse-res M17"  # "NGLL5"
            label_b = "Fine-res M00"  # "NGLL7"
        else:
            label_a = label_b = None

        self.ax.plot(tr_a_plot.times(),
                     tr_a_plot.data / tr_a_plot.data.max() + step,
                     linewidth=.75,
                     c=c_ngll5, label=label_a)

        self.ax.plot(tr_b_plot.times(),
                     tr_b_plot.data / tr_b_plot.data.max() + step,
                     linewidth=.75,
                     c=c_ngll7, label=label_b)

        self.ax.annotate(tag, xy=(0, step))

    @staticmethod
    def filter_and_evaluate(tr_a, tr_b, period, function="xcorr", 
                            tshift_max=0, **kwargs):
        """
        Calculate the value of some objective function for a given period band

        :type tr_?: obspy.core.trace.Trace
        :param tr_?: the traces with data that need to be filtered and evaluated
        :type period: float
        :param period: the period at which to filter the data
        :type function: str
        :param function: the function to use to evaluate the 'misfit'
            * xcorr: zero-lag cross correlation coefficient
            * tshift: time shift for maximum cross correlation related to a max
                shift value of 'tshift_max'
            * l1_norm: the L1 norm of the entire waveform
        :type tshift_max: float
        :param tshift_max: maximum allowable time shift in sec. for the function 
            'tshift' and 'xcorr', defaults to 0s for default function 'xcorr'
        """
        tr_a_filt = tr_a.copy()
        tr_b_filt = tr_b.copy()

        # Filter using a lowpass filter
        tr_a_filt.filter("lowpass", freq=1 / period, corners=4)
        tr_b_filt.filter("lowpass", freq=1 / period, corners=4)

        if function in ["xcorr", "tshift"]:
            cc = correlate(tr_a_filt.data, tr_b_filt.data, 
                           shift=int(tshift_max / tr_a_filt.stats.delta)
                           )
            shift, value = xcorr_max(cc)
            # If tshift requested, return the time shift value rather than the
            # maximum cross correlation coefficient
            if function == "tshift":
                value = shift * tr_a_filt.stats.delta
        elif function == "l1_norm":
            value = (np.abs(tr_a.data - tr_b.data).sum() /
                     np.sum(np.abs(tr_b.data))
                     )

        return tr_a_filt, tr_b_filt, value

    def process_single(self, code, min_period, max_period, dt, step=1.1,
                       plot=True, components=None, **kwargs):
        """
        Process a single station/component by reading in synthetic waveforms,
        filtering, evaluating objective functions and  collecting information
        into a Pandas dataframe.

        :type code: str
        :param code: station code with wildcards, e.g. NN.SSSS.?X?
        :type min_period: float
        :param min_period: minimum filter period
        :type max_period: float
        :param max_period: maximum filter period
        :type dt: float
        :param dt: step lengths to take between min and max period
        :type step: float
        :param step: for plotting, the step length between traces
        :type plot: bool
        :param plot: plot waveforms as you go
        :type components: list of str
        :param components: select which components are processed, if None,
            will check all available.
        """
        # Poplate Stream objects with all available synthetics
        st_a = Stream()
        for fid in glob(os.path.join(self.path_ngll5, f"{code}.sem?")):
            st_a += read_sem(fid, origintime=dummy_time, precision=4)

        st_b = Stream()
        for fid in glob(os.path.join(self.path_ngll7, f"{code}.sem?")):
            st_b += read_sem(fid, origintime=dummy_time, precision=4)

        # Determine which components are available
        if components is None:
            components = []
            for tr in st_a:
                components.append(tr.stats.component)

        # Resample one of the streams to the sampling rate of the other
        # ASSUMING all the sampling rates are the same within a single stream
        if st_a[0].stats.sampling_rate != st_b[0].stats.sampling_rate:
            sampling_rate = min(st_a[0].stats.sampling_rate,
                                st_b[0].stats.sampling_rate)
            for st in [st_a, st_b]:
                if st[0].stats.sampling_rate > sampling_rate:
                    st.resample(sampling_rate)

        # Incase the two traces have different lengths, trim to the shorter
        if st_a[0].stats.endtime != st_b[0].stats.endtime:
            endtime = min(st_a[0].stats.endtime, st_b[0].stats.endtime)
            for st in [st_a, st_b]:
                st.trim(dummy_time, endtime)

        # Apply filter and calculate the function evaluation
        current_step = 0
        for i, comp in enumerate(components):
            tr_a = st_a.select(component=comp)[0]
            tr_b = st_b.select(component=comp)[0]
            if plot:
                self.plot_traces(tr_a, tr_b, current_step, c_ngll7=f"C{i}",
                                 tag=f"RAW_{comp}")

            for period in np.arange(min_period, max_period, dt):
                current_step += step
                tr_a_, tr_b_, value = self.filter_and_evaluate(tr_a, tr_b,
                                                               period, **kwargs)
                # Set the stats as a series of lists in a dict
                self.stats["code"].append(self._current)
                self.stats["comp"].append(comp)
                self.stats["period"].append(period)
                self.stats["value"].append(value)

                if plot:
                    self.plot_traces(tr_a_, tr_b_, current_step, 
                                     c_ngll7=f"C{i}",
                                     tag=(f"T={period:.2f}s, "
                                          f"cc={value:.3f}")
                                     )
            current_step += step

        self.finalize_plot(**kwargs)

    def calc(self, df_fid="dfmrp.csv", **kwargs):
        """
        Main function to calculate MRP for all available stations
        
        :type df_fid: str
        :param df_fid: path and file name to save the Pandas dataframe
        """
        for code in self.stations:
            self.initiate_plot()

            net, sta, cha = code.split(".")
            self._current = f"{net}.{sta}"
            print(self._current)
            self.process_single(code, **kwargs)

        self.dfstats = pd.DataFrame(self.stats)
        self.dfstats.to_csv(df_fid, index=False)


if __name__ == "__main__":
    # Parameters
    min_period = 4
    max_period = 9.1
    dt = 1
    step = 2.25
    show = False
    save = "./figures"
    ngll5 = "./model_init"
    ngll7 = "./forest_4s"

    if save and not os.path.exists(save):
        os.makedirs(save)

    # Call Minresper
    try:
        mr = Minresper(path_ngll5=ngll5, path_ngll7=ngll7)
        mr.calc(min_period=min_period, max_period=max_period, step=step, dt=dt, 
                function="xcorr", tshift_max=0, show=show, save=save)
    except Exception as e:
        traceback.print_exc()
        import ipdb; ipdb.set_trace()

