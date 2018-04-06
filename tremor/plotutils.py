"""plotting utilities for pyfreqscan
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from utils import z2nan

mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['lines.markersize'] = 1.75


def plot_arrays(st,code_set,TEORRm,sig,show=True,save=False):
    """plot 6 streams in 3 subplot figure
    """
    f,(ax1,ax2,ax3) = plt.subplots(3,sharex=True,sharey=False,
                                                figsize=(11.69,8.27),dpi=200)

    # divy out arrays
    T,E,O,R,Rm = TEORRm
    sig2,sig3 = sig

    # create time axis
    stats = st[0].stats
    start = stats.starttime
    end = stats.endtime
    tracelength = (end-start)/3600
    t = np.linspace(0,tracelength,len(R))
    tRm = np.linspace(0,tracelength,len(Rm))

    # full waveforms
    A1 = ax1.plot(t,st[0].data,color='orange',label=st[0].get_id())
    A2 = ax1.plot(t,st[1].data,color='b',label=st[1].get_id())

    # filtered waveforms
    # B1 = ax2.plot(t,T,color='k',label="[T]remor Band 2-8Hz")
    # B2 = ax2.plot(t,E,color='r',label="[E]arthquake Band 10-20Hz")
    # B3 = ax2.plot(t,O,color='g',label="[O]cean Band .5-1Hz")
    B4 = ax2.plot(t,R,color='c',label="[R]atio")

    # amplitude ratio and median value
    # C1 = ax3.plot(t,R,color='k',label='Amplitude [R]atio (T^2/(E*O))')
    C2 = ax3.plot(tRm,Rm,color='k',label='Amplitude [R]atio [m]edian',
                                                    linewidth=1,zorder=4)
    C3 = ax3.scatter(tRm,z2nan(sig2),c='r',marker='^',zorder=5,label='2-sigma',
                                                        s=4)
    C4 = ax3.scatter(tRm,z2nan(sig3),c='orange',marker='o',zorder=6,
                                                            label='3-sigma')

    # plot vertical lines on each subplot for tremor locations
    for X,Y2 in zip(tRm,sig2):
        if np.isnan(Y2):
            continue
        else:
            for ax in [ax1,ax2,ax3]:
                D2 = ax.axvline(X,zorder=1,linestyle="solid",color='gray',
                                                        linewidth=1,alpha=0.35)

    # set axis properties
    for ax in [ax1,ax2,ax3]:
        __pretty_grids(ax)
        ax.legend(prop={"size":7.5})

    ax1.set_xlim([t.min(),t.max()])
    st_std = 5 * np.std(st[0].data)
    # ax1.set_ylim([-st_std,st_std])
    ax3.set_ylim([Rm.min(),Rm.max()])
    # ax2.set_yscale("log")

    ax1.set_title("TEROR {}".format(code_set))
    ax1.set_ylabel("velocity (m/s)")
    ax2.set_ylabel("velocity (m/s)")
    ax3.set_ylabel("dimensionless")
    ax3.set_xlabel("time (hours from UTC midnight)")
    plt.tight_layout()

    if show:
        plt.show()
    if save:
        fig_path = pathnames()['plots'] + 'tremor/'
        fig_name = code_set.format(c='?') + '.png'
        fig_out = os.path.join(fig_path,fig_name)
        plt.savefig(fig_out)
    return f


def stacked_plot(x,north_list,east_list,Rm_list,sig_list,sta_list):
    """multiple stations plotted in a stacked plot to show waveform coherence
    sts should be a stream of all components
    """
    f,(ax1,ax2,ax3) = plt.subplots(3,sharex=True,sharey=False,
                                                figsize=(11.69,8.27),dpi=200)

    # create t axis to match x
    x_Rm = np.linspace(x.min(),x.max(),len(Rm_list[0]))
    step_up = 0
    for N,E,Rm,O,S in zip(north_list,east_list,Rm_list,sig_list,sta_list):
        # full waveforms, stepping up amplitude by .5 each plot for separation
        N += step_up
        E += step_up
        step_up += 0.5
        ax1.plot(x,N)
        ax2.plot(x,E)
        # plot sigmas and Rm
        ax3.plot(x_Rm,Rm,label=S)
        ax3.scatter(x_Rm,O,marker='o')

    ax3.legend(prop={"size":7.5})
    for ax in [ax1,ax2,ax3]:
        __pretty_grids(ax)

    plt.tight_layout()
    plt.show()


    return f


def __pretty_grids(input_ax):
    """grid formatting
    """
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
                    linestyle='-',
                    linewidth='0.5',
                    color='k',
                    alpha=0.15)
    input_ax.ticklabel_format(style='sci',
                            axis='y',
                            scilimits=(0,0))
