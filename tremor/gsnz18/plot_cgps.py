import os
import glob
import matplotlib.pyplot as plt
from scipy import stats, signal
from numpy import median


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


def simple_plot():
    files_ = glob.glob("/Users/chowbr/Documents/subduction/spectral/"
                       "tremor/gsnz18/trimmed_pm50days/*_?_*")
    # files_ = glob.glob("/Users/chowbr/Documents/subduction/spectral/"
    #                    "tremor/gsnz18/trimmed_pm50days/*DNVK*")
    for fid in files_:
        sta = os.path.basename(fid).split('_')[0]
        comp = os.path.basename(fid).split('_')[1]
        with open(fid, 'r') as f:
            lines = f.readlines()
        datetime, data = [], []
        for i, line in enumerate(lines):
            i-=50
            datetime.append(i)
            data.append(float(line.strip().split(',')[1]))
            if "2017-09-08" in line.strip().split(',')[0]:
                chiapas = i
            if "2017-09-25" in line.strip().split(',')[0]:
                fiveeight = i
        f, ax = plt.subplots(1)
        pretty_grids(ax)

        # detrending stuff
        data = signal.detrend(data, type="linear")
        pre_chiapas_median = median(data[:chiapas])
        post_chiapas_median= median(data[chiapas:])
        plt.axhline(y=pre_chiapas_median, xmin=0, xmax=0.5, color='b',
                    linestyle='--', zorder=5, label="pre-chiapas median")
        plt.axhline(y=post_chiapas_median, xmin=0.5, xmax=1.0, color='b',
                    linestyle='-.', zorder=5, label="post-chiapas median")

        plt.plot(datetime, data, 'ko-')
        plt.axvline(x=chiapas, ymin=0, ymax=1, color='r', linestyle='--',
                    label="Chiapas M8.2")
        plt.axvline(x=fiveeight, ymin=0, ymax=1, color='g', linestyle='--',
                    label="2017p723941 M5.8")

        # linear regression stuff
        curve = lin_reg(datetime, data)
        pre_chiapas_curve = lin_reg(datetime[:chiapas], data[:chiapas])
        post_chiapas_curve = lin_reg(datetime[chiapas:], data[chiapas:])
        plt.plot(datetime, curve, 'r-', label="linear regression")
        # plt.plot(datetime[:chiapas], pre_chiapas_curve, 'b-',
        #          label="pre chiapas")
        # plt.plot(datetime[chiapas:], post_chiapas_curve, 'b-.',
        #          label="post chiapas")


        # misc plot stuff
        plt.legend()
        plt.xlabel("Days")
        plt.ylabel("{} (mm)".format(comp))
        plt.title("{} - displacement from intial position, DETRENDED".format(sta))
        # plt.savefig("./figures/{s}_{c}.png".format(s=sta, c=comp))
        plt.show()
        plt.close()


simple_plot()


