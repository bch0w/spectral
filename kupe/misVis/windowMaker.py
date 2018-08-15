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

def normalize_zero_to_one(array):
    """normalize an array from zero to one for easy plotting
    """
    z = (array-array.min()) / (array.max()-array.min())
    
    return z

def align_yaxis(ax1,ax2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1
    """
    ymin_a1,ymax_a1= ax1.get_ylim()
    ymin_a2,ymax_a2= ax2.get_ylim()
    
    _, y1 = ax1.transData.transform((0, (ymax_a1+ymin_a1)/2))
    _, y2 = ax2.transData.transform((0, (ymax_a2+ymin_a2)/2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    ax2.set_ylim(ymin_a2+dy, ymax_a2+dy)
    
    
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
        bounds = (-0.05,1.05)
    input_ax.set_ylim(bounds)

def setup_plot(number_of,twax=True):
    """dynamically set up plots according to number of files
    """
    # if not f:
    #     f = plt.figure(figsize=(11.69,8.27),dpi=100)

    nrows,ncols=number_of,1
    height_ratios = [1] * (number_of)
    gs = mpl.gridspec.GridSpec(nrows,ncols,height_ratios=height_ratios,hspace=0)

    # create gridspec subplots, sharex with the first axis
    axes,twaxes = [],[]
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
    axes,twaxes = setup_plot(number_of=NUMBER_OF_TRACES,twax=True)
    t = make_t_axis(st)
    complist = create_component_list(st)

    UNIT_DICT = {"DISP":"displacement [m]",
                 "VEL":"velocity [m/s]",
                 "ACC":"acceleration [m/s^2]"}
        
    for i,comp in enumerate(complist):
        # distribute data
        obs,syn = breakout_stream(st.select(component=comp))

        # plot waveforms (synthetics red, observations black)
        A1, = axes[i].plot(t,obs[0].data,'k',zorder=5,
                     label="{} (OBS)".format(obs[0].get_id()))
        A2, = axes[i].plot(t,syn[0].data,'r',zorder=5,
                     label="{} (SYN)".format(syn[0].get_id()))
        lines_for_legend = [A1,A2]

        # plot stalta data and water level
        T1, = twaxes[i].plot(t,staltas[comp],'gray',alpha=0.4,zorder=4)
        T2 = twaxes[i].axhline(y=PD["stalta_wl"],xmin=t[0],xmax=t[-1],
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
            _adj_src = normalize_zero_to_one(adj_src[comp].adjoint_source)
            T3, = twaxes[i].plot(t,_adj_src[::-1],'g',alpha=0.5,linestyle='-.',
                                label="Adjoint Source, Misfit={:.4f}".format(
                                                        adj_src[comp].misfit))
            lines_for_legend += [T3]

        except KeyError:
            pass

        # label units and twin axis
        if i == MIDDLE_TRACE:
            twaxes[i].set_ylabel("Adjoint Source (Normalized) &\n"
                                "STA/LTA, waterlevel=".format(PD["stalta_wl"]),
                                rotation=270,labelpad=20)
            comp = "{}\n{}".format(UNIT_DICT[PD['units']],comp)

        # combine legends
        labels = [l.get_label() for l in lines_for_legend]
        axes[i].legend(lines_for_legend,labels,prop={"size":9})
        
        # final adjustments
        # twaxes[i].set_yticklabels([])
        axes[i].set_ylabel(comp)
        axes[i].set_xlim([t[0],t[-1]])
        for AX in [axes[i],twaxes[i]]:
            format_axis(AX)
        align_yaxis(axes[i],twaxes[i])


            
    # figure settings
    titletext = "{s} [{b0},{b1}]".format(s=PD['station'],
                                         b0=PD['bounds'][0],
                                         b1=PD['bounds'][1])
    axes[0].set_title(titletext)
    axes[-1].set_xlabel("time [s]")

    return axes, twaxes

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
