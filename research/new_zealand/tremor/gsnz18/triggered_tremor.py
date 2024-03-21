from obspy import Stream, UTCDateTime, read
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.size'] = 15
mpl.rcParams['lines.linewidth'] = 1.
mpl.rcParams['lines.markersize'] = 1.75
mpl.rcParams['axes.linewidth'] = 2.5


def pretty_grids(input_ax):
    """make dem grids pretty
    """
    import matplotlib.ticker as ptick
    input_ax.set_axisbelow(True)
    input_ax.tick_params(which='both', direction='in', top=True, right=True)
    input_ax.minorticks_on()
    input_ax.grid(which='minor', linestyle=':', linewidth='0.5', color='k',
                  alpha=0.25)
    input_ax.grid(which='major', linestyle=':', linewidth='0.5', color='k',
                  alpha=0.25)
    # input_ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))


def setup_plot(number_of, twax=True):
    """
    Dynamically set up plots according to number_of given

    :type number_of: int
    :param number_of: number of subplots t ogenerate using gridspec
    :type twax: bool
    :param twax: whether or not twin x axes are required
    :rtype (tw)axes: matplotlib axes
    :return (tw)axes: axis objects
    """
    nrows, ncols = number_of, 1
    height_ratios = [1] * number_of
    gs = mpl.gridspec.GridSpec(nrows, ncols, height_ratios=height_ratios,
                               hspace=0)
    axes, twaxes = [], []
    for i in range(number_of):
        if i == 0:
            ax = plt.subplot(gs[i])
        else:
            ax = plt.subplot(gs[i],sharex=axes[0])
        if twax:
            twinax = ax.twinx()
            twaxes.append(twinax)
        else:
            twax = None
        pretty_grids(ax)
        axes.append(ax)

    # remove x-tick labels except for last axis
    for ax in axes[0:-1]:
        plt.setp(ax.get_xticklabels(), visible=False)

    return axes, twaxes


def origintimes():
    origintimes = [UTCDateTime("2017-09-08T23:15:41.485Z"),
                   UTCDateTime("2017-09-08T23:09:43.931Z"),
                   UTCDateTime("2017-09-08T20:43:51.334Z"),
                   UTCDateTime("2017-09-08T18:11:45.201Z"),
                   UTCDateTime("2017-09-08T18:04:58.506Z"),
                   UTCDateTime("2017-09-08T17:54:35.403Z"),
                   UTCDateTime("2017-09-08T16:56:22.540Z"),
                   UTCDateTime("2017-09-08T16:33:36.830Z"),
                   UTCDateTime("2017-09-08T15:54:16.203Z"),
                   UTCDateTime("2017-09-08T12:54:45.854Z"),
                   UTCDateTime("2017-09-08T12:07:21.307Z"),
                   UTCDateTime("2017-09-08T11:30:23.339Z"),
                   UTCDateTime("2017-09-08T11:05:55.585Z"),
                   UTCDateTime("2017-09-08T10:37:39.369Z"),
                   UTCDateTime("2017-09-08T09:56:55.879Z"),
                   UTCDateTime("2017-09-08T08:18:07.081Z"),
                   UTCDateTime("2017-09-08T07:42:47.863Z"),
                   UTCDateTime("2017-09-08T07:14:17.245Z"),
                   UTCDateTime("2017-09-08T05:34:43.632Z"),
                   UTCDateTime("2017-09-08T05:14:12.574Z"),
                   UTCDateTime("2017-09-08T04:18:34.160Z"),
                   UTCDateTime("2017-09-08T04:03:33.118Z"),
                   UTCDateTime("2017-09-08T02:53:21.334Z")
                   ]
    def convert_time_to_dec(utcdatetime):
        return utcdatetime.hour*60 + utcdatetime.minute + utcdatetime.second/60
    list_out = []
    for otime in origintimes:
        list_out.append(convert_time_to_dec(otime))
    return list_out


def preprocess(st, freqmin, freqmax):
    stp = st.copy()
    stp.trim(UTCDateTime('2017-09-08T00:00:00'),
             UTCDateTime('2017-09-09T00:00:00')
             )
    stp.decimate(2)
    stp.detrend("demean")
    stp.detrend("linear")
    stp.taper(max_percentage=0.05)
    stp.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
    # st.filter('bandpass', freqmin=5, freqmax=15)
    stp.detrend("demean")
    stp.detrend("linear")
    stp.taper(max_percentage=0.05)

    return stp


def multi_plot():
    fmin=2
    fmax=8
    pathdir = ("/seis/prj/fwi/bchow/mseeds/"
               "CHIAPAS/??.{sta}.10.?H{comp}.D.2017.251")
    stalist = ["ANWZ", "RD02", "PRHZ", "RD01", "RD03", "RD05", "RD16", "BKZ"]
    n = len(stalist)
    axes, _ = setup_plot(n+2, twax=False)
    otimes = origintimes()
    for i, sta in enumerate(stalist):
        comp = "E"
        fid_in = glob.glob(
            pathdir.format(sta=sta, comp=comp)
        )
        print(fid_in)
        if not fid_in:
            print('here')
            continue
        st = read(fid_in[0])
        print(len(st[0].data))
        st = preprocess(st, freqmin=fmin, freqmax=fmax)
        t = np.linspace(0, 1440, len(st[0].data))
        if i == 0:
            f_,ax_ = plt.subplots(1)
            plt.sca(ax_)
            plt.plot(t, st[0].data, linewidth=0.5, c='k')
        plt.sca(axes[i])
        plt.plot(t, st[0].data, linewidth=0.5, c='k')
        # plt.yticks([0], ["{}\nmax: {:.1e}m/s".format(
        #     st[0].get_id(),st[0].data.max())])
        if i != 1:
            plt.yticks([0], ["{}".format(st[0].get_id())])
        else:
            plt.yticks([0], ["{}\n(Colocated ANWZ)".format(st[0].get_id())])
        plt.axvline(x=289, ymin=0, ymax=1, linewidth=0.5, c='r')
        for otime in origintimes():
            plt.axvline(x=otime, ymin=0, ymax=1, linewidth=0.5, c='r',
                        linestyle='--')

    # surface wave band
    fid = glob.glob(pathdir.format(sta="BKZ", comp="?"))
    st = Stream()
    for f in fid:
        st += read(f)
    st = preprocess(st, freqmin=1/200, freqmax=1/60)
    for ax, comp in zip([-1,-2], ["R", "T"]):
        st_s = st.select(component=comp)
        plt.sca(axes[ax])
        plt.plot(t, st_s[0].data, linewidth=0.5, c='k')
        plt.yticks([0], ["{}\n(60-200s)".format(st_s[0].get_id())])
        plt.axvline(x=289, ymin=0, ymax=1, linewidth=0.5, c='r')
        for otime in origintimes():
            plt.axvline(x=otime, ymin=0, ymax=1, linewidth=0.5, c='r',
                        linestyle='--')

    plt.sca(axes[0])
    # plt.title("Southern Hawke's Bay [{}-{}Hz Bandpass]".format(fmin,fmax))
    plt.sca(axes[-1])
    plt.xlabel('2017.251 (minutes since 00:00:00.0Z)')
    plt.xlim([t.min(), t.max()])
    plt.show()


if __name__ == "__main__":
    multi_plot()

    # gps2dist_azimuth(-40.6796,176.2462,15.068,-93.715)
