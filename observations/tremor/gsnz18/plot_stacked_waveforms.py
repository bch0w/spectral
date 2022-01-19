import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
import matplotlib as mpl

from obspy import read


mpl.rcParams['font.size'] = 17.5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['lines.markersize'] = 1.75
mpl.rcParams['axes.linewidth'] = 1.25


def create_min_max(tr,pixel_length=10):
    """
    Creates new data using a min/max approach that calculated the minimum and
    maximum values of each "pixel". Much faster plotting of large datasets.

    !!! base code copied and modified
    !!! from obspy.imaging.waveform.__plot_min_max
    """
    # Some variables to help calculate the values.
    starttime = date2num(tr.stats.starttime.datetime)
    endtime = date2num(tr.stats.endtime.datetime)
    # The same trace will always have the same sampling_rate.
    sampling_rate = tr.stats.sampling_rate
    # width of x axis in seconds
    x_width = endtime - starttime
    # number of samples that get represented by one min-max pair
    # width = 800 # guessing
    # pixel_length = int(
    #     np.ceil((x_width * sampling_rate + 1) / width))
    # Loop over all the traces. Do not merge them as there are many samples
    # and therefore merging would be slow.
    trace_length = len(tr.data)
    pixel_count = int(trace_length // pixel_length)
    remaining_samples = int(trace_length % pixel_length)
    remaining_seconds = remaining_samples / sampling_rate
    # Reference to new data array which does not copy data but can be
    # reshaped.
    if remaining_samples:
        data = tr.data[:-remaining_samples]
    else:
        data = tr.data
    data = data.reshape(pixel_count, pixel_length)
    min_ = data.min(axis=1) * tr.stats.calib
    max_ = data.max(axis=1) * tr.stats.calib
    # Calculate extreme_values and put them into new array.
    if remaining_samples:
        extreme_values = np.empty((pixel_count + 1, 2), dtype=np.float)
        extreme_values[:-1, 0] = min_
        extreme_values[:-1, 1] = max_
        extreme_values[-1, 0] = \
            tr.data[-remaining_samples:].min() * tr.stats.calib
        extreme_values[-1, 1] = \
            tr.data[-remaining_samples:].max() * tr.stats.calib
    else:
        extreme_values = np.empty((pixel_count, 2), dtype=np.float)
        extreme_values[:, 0] = min_
        extreme_values[:, 1] = max_
    # Finally plot the data.
    start = date2num(tr.stats.starttime.datetime)
    end = date2num(tr.stats.endtime.datetime)
    if remaining_samples:
        # the last minmax pair is inconsistent regarding x-spacing
        x_values = np.linspace(start, end - remaining_seconds,
                               num=extreme_values.shape[0] - 1)
        x_values = np.concatenate([x_values, [end]])
    else:
        x_values = np.linspace(start, end, num=extreme_values.shape[0])
    x_values = np.repeat(x_values, 2)
    y_values = extreme_values.flatten()

    return x_values, y_values


def clip_data(y):
    sigma = np.std(y)
    y[y > 6 * sigma] = np.nan
    y[y < -6 * sigma] = np.nan

    return y

def time_to_decimal(datetime):
    return datetime.julday + (datetime.hour + datetime.minute/60)/24


def pretty_grids(input_ax, scitick=False):
    """make dem grids pretty
    """
    import matplotlib.ticker as ptick
    input_ax.set_axisbelow(True)
    input_ax.tick_params(which='both',
                         direction='in',
                         top=True,
                         right=True)
    input_ax.minorticks_on()
    input_ax.grid(which='minor',
                    linestyle=':',
                    linewidth='0.5',
                    color='k',
                    alpha=0.25)
    input_ax.grid(which='major',
                    linestyle=':',
                    linewidth='0.5',
                    color='k',
                    alpha=0.25)
    if scitick:
        input_ax.ticklabel_format(style='sci',
                                axis='y',
                                scilimits=(0,0))


files_ = glob.glob('*?HE*')
files_ += glob.glob('../filtered_3-60s/*?HE*')
f, ax = plt.subplots(1)
pretty_grids(ax)
chiapas = 251.2
step = 0
steps, ids = [],[]
for i, fid in enumerate(files_):
    print(fid)
    st = read(fid)
    _, y = create_min_max(st[0])
    x = np.linspace(time_to_decimal(st[0].stats.starttime),
                    time_to_decimal(st[0].stats.endtime),
                    len(y)
                    )
    print(len(y))
    if i != len(files_)-1:
        y /= y.max()
        y += step

        plt.plot(x, y, linewidth=0.5, c='k', zorder=10)

    else:
        y /= y.max()
        y += step
        plt.plot(x, y, linewidth=0.5, zorder=10, c='r')

    # plt.scatter(x=x[np.argmax(y)], y=y.max(), zorder=1000, c='r')
    # plt.scatter(x=x[np.argmax(y)], y=step, zorder=1000, c='r')

    ids.append("{}\npeak: {:.1e} m/s".format(st[0].get_id(), st[0].data.max()))
    steps.append(step)
    step += 0.4

plt.title('Waveforms [filtered 2-8Hz except RD16, filtered 3-60s]')
plt.yticks(steps, ids, fontsize=12)
plt.axvline(chiapas, 0, 1, c='r', zorder=2, linewidth=1.5)
plt.xlabel('Julian Days (2017)')
plt.xlim([x.min(),x.max()])
# plt.setp(ax.get_yticklabels(), visible=False)
plt.show()
