import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from scipy import stats, signal

import matplotlib as mpl
mpl.rcParams['font.size'] = 16
mpl.rcParams['lines.linewidth'] = 1.25
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


def lin_reg(x_in, y_in):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_in, y_in)
    y_out = lambda m, x, b: m*x + b
    data_out = []
    for x in x_in:
        data_out.append(y_out(m=slope, x=x, b=intercept))

    return data_out


def simple_plot(detrend=False):
    f, ax = plt.subplots(1)
    pretty_grids(ax)
    for sta in ["OROA", "MNHR", "PORA", "AKTO", "PAWA"]:
        files_ = []
        for comp in ["e"]:#, "e"]:
            files_ += glob.glob("/seis/prj/fwi/bchow/spectral/tremor/gsnz18/cgps/"
                               "trimmed_pm60days/{sta}_{c}_*".format(
                sta=sta, c=comp))
        # files_ = glob.glob("/Users/chowbr/Documents/subduction/spectral/"
        #                    "tremor/gsnz18/trimmed_pm50days/*DNVK*")
        for fid in files_:
            sta = os.path.basename(fid).split('_')[0]
            comp = os.path.basename(fid).split('_')[1]
            with open(fid, 'r') as f:
                lines = f.readlines()
            datetime, data = [], []
            for i, line in enumerate(lines):
                datetime.append(UTCDateTime(line.strip().split(',')[0]))
                data.append(float(line.strip().split(',')[1]))
                if "2017-09-08" in line.strip().split(',')[0]:
                    chiapas = i
                if "2017-09-25" in line.strip().split(',')[0]:
                    fiveeight = i

            # detrending stuff
            if detrend:
                data = signal.detrend(data, type="linear")
                detrend = "DETRENDED"
            else:
                detrend = ""
            # pre_chiapas_median = median(data[:chiapas])
            # post_chiapas_median= median(data[chiapas:])
            # plt.axhline(y=pre_chiapas_median, xmin=0, xmax=0.5, color='b',
            #             linestyle='--', zorder=5, label="pre-chiapas median")
            # plt.axhline(y=post_chiapas_median, xmin=0.5, xmax=1.0, color='b',
            #             linestyle='-.', zorder=5, label="post-chiapas median")
            t = np.linspace(datetime[0].julday,datetime[-1].julday,len(data))
            plt.plot(t, data, 'o-', markersize=5, label=sta,
                     markerfacecolor='w', markeredgewidth=1.2, zorder=5)

    plt.axvline(x=251, ymin=0, ymax=1, color='k', linestyle='--', linewidth=2,
                label="Chiapas M8.2", zorder=2)
    #plt.axvline(x=fiveeight, ymin=0, ymax=1, color='g', linestyle='--',
    #            label="2017p723941 M5.8")

    # linear regression stuff
    curve = lin_reg(t, data)
    # pre_chiapas_curve = lin_reg(datetime[:chiapas], data[:chiapas])
    # post_chiapas_curve = lin_reg(datetime[chiapas:], data[chiapas:])
    # plt.plot(t, curve, 'r-', label="linear regression")
    # plt.plot(datetime[:chiapas], pre_chiapas_curve, 'b-',
    #          label="pre chiapas")
    # plt.plot(datetime[chiapas:], post_chiapas_curve, 'b-.',
    #          label="post chiapas")


    # misc plot stuff
    plt.legend(loc="best")
    plt.xlabel("Julian Days 2017")
    plt.ylabel("{} (mm)".format(comp))
    plt.title("cGPS displacement from intial position")
    # plt.savefig("./figures/{s}_{c}.png".format(s=sta, c=comp))
    plt.show()
    plt.close()


if __name__ == "__main__":
    simple_plot()


