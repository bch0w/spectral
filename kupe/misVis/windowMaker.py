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
        
def window_maker(st,windows,*args,**kwargs):
    """plot streams and windows. assumes you have N observation traces and 
    N synthetic traces for a 2N length stream object
    """
    # function setup
    NUMBER_OF_TRACES = len(st)//2
    MIDDLE_TRACE = NUMBER_OF_TRACES//2
    axes,_ = setup_plot(number_of=NUMBER_OF_TRACES,twax=False)
    t = make_t_axis(st)
    complist = create_component_list(st)
    
    for i,comp in enumerate(complist):
        # distribute data
        obs,syn = breakout_stream(st.select(component=comp))
        windowlist = windows[comp]
        
        # plot waveforms (synthetics red, observations black)
        axes[i].plot(t,obs[0].data,'k',
                     label="{} (OBS)".format(obs[0].get_id()),
                     zorder=5)
        axes[i].plot(t,syn[0].data,'r',
                     label="{} (SYN)".format(syn[0].get_id()),
                     zorder=5)
        
        # plot windows as transparent boxes
        ymin,ymax = axes[i].get_ylim()
        for window in windowlist:
            xwindow = np.arange(window.left,window.right,1)
            twindow = t[xwindow]
            axes[i].fill_between(twindow,ymin,ymax,
                                 facecolor='orange',
                                 edgecolor='k',
                                 linewidth=0.5,
                                 alpha=0.25)
                                 
            # annotate boxes with information from window
            winanno = "maxCC:{mcc:.4f}\nccShift:{ccs}\ndlnA:{dln:.4f}".format(
                                            mcc=window.max_cc_value,
                                            ccs=window.cc_shift,
                                            dln=window.dlnA)
            axes[i].annotate(s=winanno,xy=(twindow[10],ymax*0.5),
                                                    zorder=4,fontsize=7)
            
        # final plot adjustments
        axes[i].set_xlim([t[0],t[-1]])
        axes[i].legend(prop={"size":9})
        if i == MIDDLE_TRACE:
            comp = "displacement [m]\n{}".format(comp)
        axes[i].set_ylabel(comp)
    
    # figure settings
    PD = kwargs['PD']
    
    titletext = "{s} [{b0},{b1}]".format(s=PD['station_name'],
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
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    