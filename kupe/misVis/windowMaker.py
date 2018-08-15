"""misfit visualization tool to be called through adjointBuilder
produces waveform plots for a stream object and plots misfit windows
outputted by pyflex
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from adjointBuilder import breakout_stream

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.linewidth'] = 2

def pretty_grids(input_ax):
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
    input_ax.ticklabel_format(style='sci',
                            axis='y',
                            scilimits=(0,0))

def format_axis(input_ax):
    """sit the tick marks away from the plot edges to prevent overlapping when
    multiple subplots are stacked atop one another, and for general gooood looks
    will check if the plot is two sided (e.g. waveforms) or only positive
    (e.g. STA/LTA)
    """
    ymin,ymax= input_ax.get_ylim()
    maxvalue = max([abs(_) for _ in input_ax.get_ylim()])
    percentover = maxvalue*0.125
    if abs(round(ymin/ymax)) == 1:
        bounds = (-1*(maxvalue+percentover),(maxvalue+percentover))
    else:# elif abs(round(ymin/ymax)) == 0:
        bounds = (-0.05,(maxvalue+percentover))
    input_ax.set_ylim(bounds)

def setup_plot(number_of,num_twax=0):
    """dynamically set up plots according to number of files
    twax sets up twin x axes for separate plotting, integer value for number of
    extra twin x axes to set. if twax=0 twax returns an empty list, twaxes
    """
    nrows,ncols=number_of,1
    height_ratios = [1] * (number_of)
    gs = mpl.gridspec.GridSpec(nrows,ncols,height_ratios=height_ratios,hspace=0)

    # create gridspec subplots, sharex with the first axis
    axes = []
    twaxes = [[] for _ in range(num_twax)]
    for i in range(number_of):
        if i == 0:
            ax = plt.subplot(gs[i])
        else:
            ax = plt.subplot(gs[i],sharex=axes[0])
        # dynamically generate twinx axes as a list of lists
        twinax = ax.twinx()
        for t in range(num_twax):
            twaxes[t].append(twinax)

        pretty_grids(ax)
        axes.append(ax)

    # remove x-tick labels except for last axis
    for ax in axes[0:-1]:
        plt.setp(ax.get_xticklabels(),visible=False)

    return axes,twaxes

def make_t_axis(st):
    """return a time axis for plotting
    """
    stats = st[0].stats
    t = np.linspace(0,stats.endtime-stats.starttime,len(st[0].data))
    return t

def create_component_list(st):
    """figure out if NEZ or RTZ components for proper list iteration
    """
    complist = []
    for tr in st:
        complist.append(tr.get_id()[-1])
    complist = list(set(complist))
    complist.sort()

    return complist

def window_maker(st,windows,staltas,adj_src,*args,**kwargs):
    """plot streams and windows. assumes you have N observation traces and
    N synthetic traces for a 2N length stream object
    """
    PD = kwargs['PD']

    # function setup
    NUMBER_OF_TRACES = len(st)//2
    MIDDLE_TRACE = NUMBER_OF_TRACES//2
    axes,twaxes = setup_plot(number_of=NUMBER_OF_TRACES,num_twax=2)
    t = make_t_axis(st)
    complist = create_component_list(st)

    UNIT_DICT = {"DISP":"displacement [m]",
                 "VEL":"velocity [m/s]",
                 "ACC":"acceleration [m/s^2]"}

    # distribute twinaxes for STALTA and adjoint source plotting, 
    # put adjoint source axes infront of STALTA
    ax_stalta = twaxes[0]
    ax_adjsrc = twaxes[1]
    for _ax_stalta,_ax_adjsrc in zip(ax_adjsrc,ax_stalta):
        _ax_adjsrc.set_zorder(_ax_stalta.get_zorder()+1)
        # _ax_adjsrc.set_visible(False)

    for i,comp in enumerate(complist):
        # distribute data
        obs,syn = breakout_stream(st.select(component=comp))

        # plot waveforms (synthetics red, observations black)
        A1, = axes[i].plot(t,obs[0].data,'k',
                     label="{} (OBS)".format(obs[0].get_id()),
                     zorder=5)
        A2, = axes[i].plot(t,syn[0].data,'r',
                     label="{} (SYN)".format(syn[0].get_id()),
                     zorder=5)

        # plot stalta data and water level
        T1, = ax_stalta[i].plot(t,staltas[comp],'gray',alpha=0.4,zorder=4,
                        label="STA/LTA (SYN), WL={}".format(PD["stalta_wl"]))
        T2 = ax_stalta[i].axhline(y=PD["stalta_wl"],xmin=t[0],xmax=t[-1],
                      alpha=0.2,zorder=3,linewidth=1.5,c='k',linestyle='--')
        
        # plot windows (if available) as semi-transparent boxes
        ymin,ymax = axes[i].get_ylim()
        window_anno_template = "maxCC:{mcc:.4f}\nccShift:{ccs}\ndlnA:{dln:.4f}"
        try:
            for window in windows[comp]:
                xwindow = np.arange(window.left,window.right,1)
                twindow = t[xwindow]
                F1 = axes[i].fill_between(twindow,ymin,ymax,
                                     facecolor='orange',
                                     edgecolor='k',
                                     linewidth=0.5,
                                     alpha=0.25)

                # annotate boxes with information from window
                winanno = window_anno_template.format(mcc=window.max_cc_value,
                                                      ccs=window.cc_shift,
                                                      dln=window.dlnA)
                axes[i].annotate(s=winanno,xy=(twindow[10],ymax*0.5),
                                                        zorder=4,fontsize=7)
                                                        
            # plot adjoint source if available - if there are windows, then 
            # there is also an adjoint source to plot
            T3, = ax_adjsrc[i].plot(t,adj_src[comp].adjoint_source,
                                    'g',alpha=0.75,
                                    label="Adjoint Source, Misfit={}".format(
                                                        adj_src[comp].misfit))

        except KeyError:
            pass

        # label units and twin axis
        if i == MIDDLE_TRACE:
            ax_adjsrc[i].set_ylabel("Adjoint Source & STA/LTA",rotation=90)
            comp = "{}\n{}".format(UNIT_DICT[PD['units']],comp)

            # combine legends
            lines = [A1,A2,T1,T3]
            labels = [l.get_label() for l in lines]
            axes[i].legend(lines,labels,prop={"size":9})

        axes[i].set_ylabel(comp)
        axes[i].set_xlim([t[0],t[-1]])
        for AX in [axes[i],ax_stalta[i]]:
            format_axis(AX)


    # figure settings
    titletext = "{s} [{b0},{b1}]".format(s=PD['station'],
                                         b0=PD['bounds'][0],
                                         b1=PD['bounds'][1])
    axes[0].set_title(titletext)
    axes[-1].set_xlabel("time [s]")

    return axes

def _test_window_maker():
    """check that window maker works proper by using internal test data
    """
    boundsdict = {"station_name":"TEST","bounds":[(6,30)]}
    streampath = pathnames()['data'] + 'WINDOWTESTING/testmseed.pickle'
    windowpath = pathnames()['data'] + 'WINDOWTESTING/testwindows.npz'
    st = read(streampath)
    windows = np.load(windowpath)
    axes = window_maker(st,windows,PD=boundsdict)


if __name__ == "__main__":
    _test_window_maker()
