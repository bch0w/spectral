"""misfit visualization tool to be called through adjointBuilder
produces waveform plots for a stream object as well
"""
import matplotlib as mpl
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.linewidth'] = 1.5

# =============================== HELPER FUNCTIONS ============================
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
                            
def setup_plot(number_of_files,twax=True):
    """dynamically set up plots according to number of files
    """
    f = plt.figure(figsize=(11.69,8.27),dpi=75)
    nrows,ncols=number_of_files,1
    height_ratios = [1] * (number_of_files)
    gs = gridspec.GridSpec(nrows,ncols,height_ratios=height_ratios,hspace=0)

    # create gridspec subplots, sharex with the first axis
    axes,twaxes = [],[]
    for i in range(number_of_files):
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

    # remove x tick labels except for last axis
    for ax in axes[0:-1]:
        plt.setp(ax.get_xticklabels(),visible=False)
    
    return f,axes,twaxes
    
def windowMaker(st,windows):
    """plot streams and windows
    """
    
    
    
    
    