"""plotting utilities for pyfreqscan
"""
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec

# internal packages
sys.path.append("../modules")
from getdata import pathnames

from utils import __z2nan

mpl.rcParams['font.size'] = 10
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['lines.markersize'] = 1.75
mpl.rcParams['axes.linewidth'] = 1.25


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
    C3 = ax3.scatter(tRm,__z2nan(sig2),c='r',marker='^',zorder=5,
                                                            label='2-sigma',s=4)
    C4 = ax3.scatter(tRm,__z2nan(sig3),c='orange',marker='o',zorder=6,
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

def stacked_plot(code,x,N,E,**options):
    """multiple stations plotted in a stacked plot to show waveform coherence
    sts should be a stream of all components
    """
    net,sta,loc,cha,d,year,jday = code.split('.')

    f,(ax1,ax2,ax3) = plt.subplots(3,sharex=True,sharey=False,
                                        figsize=(11.69,8.27),dpi=75)
    ax3a = ax3.twinx()

    # create t axis to match x
    TEORRm_ = options['TEORRm_list'][0]
    x_Rm = np.linspace(x.min(),x.max(),len(TEORRm_['Rm']))
    x_Rs = np.linspace(x.min(),x.max(),len(TEORRm_['Rs']))
    step_up = 0
    n_ano,e_ano = options['ano_list']
    for N_,E_,T_,TEROR_,O_,Nano_,Eano_ in zip(N,E,
                                            options['tremor_list'],
                                            options['TEORRm_list'],
                                            options['sig_list'],
                                            n_ano,
                                            e_ano):
        Rm_ = TEROR_['Rm']
        Rs_ = TEROR_['Rs']

        N_ += step_up
        E_ += step_up
        T_ += step_up
        ax1.set_title('[Tremor Detection] {}\n \
                    North (2-8Hz Bandpass, 0-1 norm, 0.5 cutoff)'.format(jday))
        ax1.plot(x,N_)
        ax1.plot(x,T_,'k')
        ax2.set_title('East')
        ax2.plot(x,E_)
        ax2.plot(x,T_,'k')
        ax3.set_title('Ratio median detection threshold and 2-sigma detection')
        ax3.plot(x_Rm,Rm_,'|-',markersize=1.25)
        # plot Rs list on twin ax3
        ax3a.plot(x_Rs,Rs_,'|-',markersize=1.25)
        ax3.scatter(x_Rm,O_,marker='o',s=4,zorder=10)
        ax3.axhline(y=options['horizontal_line'],xmin=x.min(),xmax=x.max(),c='k')
        ax3.set_xlabel('NZ local time (hours)')
        for ax,ano in zip([ax1,ax2],[Nano_,Eano_]):
            ax.annotate(ano,xy=(x.min(),step_up),
                                        xytext=(x.min(),step_up),fontsize=6)
        step_up += 0.25

    # ax3.legend(prop={"size":7.5})
    for ax in [ax1,ax2,ax3,ax3a]:
        __pretty_grids(ax)
        # ax.set_xlim([18,30]) if night else ax.set_xlim([x.min(),x.max()])
        ax.set_xlim([x.min(),x.max()])

    plt.tight_layout()
    if options['save']:
        fig_path = pathnames()['plots'] + 'tremor/'
        fig_name = code_set.format(c='?') + '.png'
        fig_out = os.path.join(fig_path,fig_name)
        plt.savefig(fig_out,dpi=200)

    if options['show']:
        plt.show()

    plt.close()

def gridspec_plot(code,x,N,E,**options):
    """based on stacked_plot, but with gridspec changing size of subplots:
    multiple stations plotted in a stacked plot to show waveform coherence
    sts should be a stream of all components
    
    gridspec_plot(code=code_set,x=t,N=y_N_list,E=y_E_list,
                 TEORRm_list=TEROR_list,
                 sig_list=sig_list,
                 ano_list=ano_NE_list,
                 tremor_list=tremor_list,
                 night=night,
                 show=True,
                 save=False)
    """
    net,sta,loc,cha,d,year,jday = code.split('.')

    # initiate figure and axes objects
    f = plt.figure(figsize=(11.69,8.27),dpi=75)
    gs = gridspec.GridSpec(5,1,height_ratios=[2,2,1,1,1],hspace=0.2)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1],sharex=ax1)
    ax3a = plt.subplot(gs[2],sharex=ax1)
    ax3b = plt.subplot(gs[3],sharex=ax1)
    ax3c = plt.subplot(gs[4],sharex=ax1)
    for ax in [ax1,ax2,ax3a,ax3b]:
        plt.setp(ax.get_xticklabels(), visible=False)


    # create t axis to match x
    TEORRm_ = options['TEORRm_list'][0]
    x_Rm = np.linspace(x.min(),x.max(),len(TEORRm_['Rm']))
    x_Rs = np.linspace(x.min(),x.max(),len(TEORRm_['Rs']))
    x_Rh = np.linspace(x.min(),x.max(),len(TEORRm_['Rh']))

    step_up = 0
    n_ano,e_ano = options['ano_list']
    for N_,E_,T_,TEROR_,TREMOR_,Nano_,Eano_ in zip(N,E,
                                            options['tremor_list'],
                                            options['TEORRm_list'],
                                            options['sig_list'],
                                            n_ano,
                                            e_ano):
        # amplitude ratios
        Rs_ = TEROR_['Rs']
        Rm_ = TEROR_['Rm']
        Rh_ = TEROR_['Rh']
        # scatterpoints for values about two-sigma
        Rs_sigma_ = TREMOR_['Rs']
        Rm_sigma_ = TREMOR_['Rm']
        Rh_sigma_ = TREMOR_['Rh']
        # shifting waveform plots up by increments
        N_ += step_up
        E_ += step_up
        T_ += step_up

        # waveforms
        ax1.plot(x,N_)
        ax1.plot(x,T_,'k')
        ax2.plot(x,E_)
        ax2.plot(x,T_,'k')
        # ratios
        ax3a.plot(x_Rs,Rs_,'|-',markersize=0.25)
        ax3b.plot(x_Rm,Rm_,'|-',markersize=1.25)
        ax3c.plot(x_Rh,Rh_,'|-',markersize=1.25)
        # two-sigmas
        ax3a.scatter(x_Rs,Rs_sigma_,marker='o',s=1.5,zorder=10)
        ax3b.scatter(x_Rm,Rm_sigma_,marker='o',s=4,zorder=10)
        ax3c.scatter(x_Rh,Rh_sigma_,marker='o',s=4,zorder=10)
        # annotations
        for ax,ano in zip([ax1,ax2],[Nano_,Eano_]):
            ax.annotate(ano,xy=(x.min(),step_up),
                        xytext=(x.min(),step_up),
                        fontsize=8,
                        zorder=11)
        step_up += 0.5

    # labelliung
    title = ("[Tremor Detection]\n"
         "Jday: {} Bandpass: [2-5Hz] Norm: 0-1 Cutoff: 0.5".format(jday))
    ax1.set_title(title)
    ax1.set_ylabel('N. velocity (-1/1 norm)')
    ax2.set_ylabel('E. velocity (-1/1 norm)')
    ax3a.set_ylabel('5-sec. ratio (R$_s$)')
    ax3b.set_ylabel('5-min. ratio (R$_m$)')
    ax3c.set_ylabel('1-hour ratio (R$_h$)')
    ax3c.set_xlabel('NZ local time (hours)')

    # axis options
    # ax3.legend(prop={"size":7.5})
    for ax in [ax1,ax2,ax3a,ax3b,ax3c]:
        __pretty_grids(ax)
        # ax.set_xlim([18,30]) if night else ax.set_xlim([x.min(),x.max()])
        ax.set_xlim([x.min(),x.max()])
    # plt.tight_layout()

    if options['save']:
        fig_path = pathnames()['plots'] + 'tremor/'
        fig_name = code.format(c='?') + '.png'
        fig_out = os.path.join(fig_path,fig_name)
        plt.savefig(fig_out,dpi=200)

    if options['show']:
        plt.show()

    plt.close()
    
def envelope_plots(code,x,N,E,**options):
    """based on stacked_plot, but with gridspec changing size of subplots:
    multiple stations plotted in a stacked plot to show waveform coherence
    sts should be a stream of all components
    """
    net,sta,loc,cha,d,year,jday = code.split('.')

    # initiate figure and axes objects
    f = plt.figure(figsize=(11.69,8.27),dpi=75)
    gs = gridspec.GridSpec(6,1,height_ratios=[2,2,2,1,1,1],hspace=0.2)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1],sharex=ax1)
    ax2a = plt.subplot(gs[2],sharex=ax1)
    ax3a = plt.subplot(gs[3],sharex=ax1)
    ax3b = plt.subplot(gs[4],sharex=ax1)
    ax3c = plt.subplot(gs[5],sharex=ax1)
    for ax in [ax1,ax2,ax3a,ax3b]:
        plt.setp(ax.get_xticklabels(), visible=False)


    # create t axis to match x
    TEORRm_ = options['TEORRm_list'][0]
    x_Rm = np.linspace(x.min(),x.max(),len(TEORRm_['Rm']))
    x_Rs = np.linspace(x.min(),x.max(),len(TEORRm_['Rs']))
    x_Rh = np.linspace(x.min(),x.max(),len(TEORRm_['Rh']))
    x_Env = np.linspace(x.min(),x.max(),len(options['envelopes'][0]))
    x_sur = np.linspace(x.min(),x.max(),len(options['surface']))

    step_up = 0
    n_ano,e_ano = options['ano_list']
    for N_,E_,Env_,TEROR_,TREMOR_,Nano_,Eano_ in zip(N,E,
                                            options['envelopes'],
                                            options['TEORRm_list'],
                                            options['sig_list'],
                                            n_ano,
                                            e_ano):
        # amplitude ratios
        Rs_ = TEROR_['Rs']
        Rm_ = TEROR_['Rm']
        Rh_ = TEROR_['Rh']
        # scatterpoints for values about two-sigma
        Rs_sigma_ = TREMOR_['Rs']
        Rm_sigma_ = TREMOR_['Rm']
        Rh_sigma_ = TREMOR_['Rh']
        # shifting waveform plots up by increments
        N_ += step_up
        # Env_ += step_up

        # waveforms
        ax1.plot(x,N_)
        ax2.plot(x_Env,Env_)
        # ratios
        ax3a.plot(x_Rs,Rs_,'|-',markersize=0.25)
        ax3b.plot(x_Rm,Rm_,'|-',markersize=1.25)
        ax3c.plot(x_Rh,Rh_,'|-',markersize=1.25)
        # two-sigmas
        ax3a.scatter(x_Rs,Rs_sigma_,marker='o',s=1.5,zorder=10)
        ax3b.scatter(x_Rm,Rm_sigma_,marker='o',s=4,zorder=10)
        ax3c.scatter(x_Rh,Rh_sigma_,marker='o',s=4,zorder=10)
        # annotations
        for ax,ano in zip([ax1,ax2],[Nano_,Eano_]):
            ax.annotate(ano,xy=(x.min(),step_up),
                        xytext=(x.min(),step_up),
                        fontsize=8,
                        zorder=11)
        step_up += 0.5
    ax2a.plot(x_sur,options['surface'],'k')


    # labelliung
    title = ("[Tremor Detection]\n"
         "Jday: {} Bandpass: [2-5Hz] Norm: 0-1 Cutoff: 0.5".format(jday))
    ax1.set_title(title)
    ax1.set_ylabel('N. velocity (-1/1 norm)')
    ax2.set_ylabel('Waveform Envelopes')
    ax2a.set_ylabel("N. Velocity (m/s) [6-30s]")
    ax3a.set_ylabel('5-sec. ratio (R$_s$)')
    ax3b.set_ylabel('5-min. ratio (R$_m$)')
    ax3c.set_ylabel('1-hour ratio (R$_h$)')
    ax3c.set_xlabel('NZ local time (hours)')

    # axis options
    # ax3.legend(prop={"size":7.5})
    for ax in [ax1,ax2,ax2a,ax3a,ax3b,ax3c]:
        __pretty_grids(ax)
        # ax.set_xlim([18,30]) if night else ax.set_xlim([x.min(),x.max()])
        ax.set_xlim([x.min(),x.max()])
    # plt.tight_layout()

    if options['save']:
        fig_path = pathnames()['plots'] + 'tremor/'
        fig_name = code.format(c='?') + '.png'
        fig_out = os.path.join(fig_path,fig_name)
        plt.savefig(fig_out,dpi=200)

    if options['show']:
        plt.show()

    plt.close()

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
                    linestyle='--',
                    linewidth='0.55',
                    color='k',
                    alpha=0.25)
    input_ax.ticklabel_format(style='sci',
                            axis='y',
                            scilimits=(0,0))
