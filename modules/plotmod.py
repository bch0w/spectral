"""module file for plotting functions
"""
def pretty_grids(input_ax):
    """make dem grids pretty
    """
    input_ax.set_axisbelow(True)
    input_ax.tick_params(which='both',direction='in',top=True,right=True)
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


def align_yaxis(ax1,v1,ax2,v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1
    """
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)
