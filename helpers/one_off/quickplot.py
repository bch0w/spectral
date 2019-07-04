"""plotting specfem outputs to check the differences
"""
import os
import glob
import matplotlib.pyplot as plt

def linespecs():
    import matplotlib as mpl
    mpl.rcParams['font.size'] = 12
    mpl.rcParams['lines.linewidth'] = 1.5
    mpl.rcParams['lines.markersize'] = 1.75
    mpl.rcParams['axes.linewidth'] = 2.0


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


def using_obspy():
    """
    if signal processing is required
    """
    from obspy import read

    base = ('/Users/chowbr/Documents/subduction/tomo/mesh_test/')
    folder_1 = "srtm30p_downsamp_ngll5"
    folder_2 = "etopo1_carl"

    output_1 = os.path.join(base,folder_1,'*.mseed')
    semd1 = glob.glob(output_1)
    semd1.sort()

    semd2 = []
    for fid1 in semd1:
        check1 = os.path.basename(fid1)
        output_2 = os.path.join(base, folder_2, check1)
        if os.path.exists(output_2):
            semd2.append(output_2)
        else:
            semd1.remove(fid1)
   
    linespecs()
    for ms1, ms2 in zip(semd1, semd2):
        print(os.path.basename(ms1), os.path.basename(ms2))
        st1 = read(ms1)
        st2 = read(ms2)

        f,ax = plt.subplots(1)
        for st, c, f in zip([st1, st2], ["r", "k"], [folder_1, folder_2]):
            st.taper(max_percentage=0.05)
            st.filter('bandpass', freqmin=1/30, freqmax=1/10)
            plt.plot(st[0].data, color=c, label=f)
        
        plt.title("{} 4km Coarse (red) vs. 1km Fine (black); 10-30s".format(st[0].get_id()))
        plt.xlabel("Samples")
        plt.ylabel("Amplitude")
        pretty_grids(ax)
        plt.show()

def single_file():
    """
    Compare single files in a shared folder
    """
    base = ('/Users/chowbr/Documents/subduction/tomo/mesh_test/2018p130600')
    output_1_name = 'mesh_original_xxrd01.semv'
    output_2_name = 'mesh_downsamp_xxrd01.semv'

    output_1 = os.path.join(base,output_1_name)
    output_2 = os.path.join(base,output_2_name)

    for O,color in zip([output_1, output_2],["r", "k"]):
        lines = open(O,'r').readlines()
        amps, times = [],[]
        for line in lines:
            time, amp = line.strip().split()
            amps.append(float(amp))
            times.append(float(time))
        plt.plot(times,amps,color,label=O)

    plt.legend()
    plt.title('{} (black) vs {} (red)'.format(output_1_name,output_2_name))
    plt.show()

def folder_compare():
    """
    Compare all files in a folder, not a smart function
    """
    base = ('/Users/chowbr/Documents/subduction/tomo/mesh_test/')
    folder_1 = "srtm30p_downsamp"
    folder_2 = "etopo1_carl"

    output_1 = os.path.join(base,folder_1,'*.semv')
    semd1 = glob.glob(output_1)
    semd1.sort()

    output_2 = os.path.join(base,folder_2,'*.semv')
    semd2 = glob.glob(output_2)
    semd2.sort()

    for Y,D in zip(semd1,semd2):
        y_amps,y_times,d_amps,d_times = [],[],[],[]
        y_fid = open(Y,'r')
        y_fid = y_fid.readlines()
        d_fid = open(D,'r')
        d_fid = d_fid.readlines()
        for YY,DD in zip(y_fid,d_fid):
            yT,yA = YY.strip().split()
            dT,dA = DD.strip().split()
            y_amps.append(float(yA))
            y_times.append(float(yT))
            d_amps.append(float(dA))
            d_times.append(float(dT))
        plt.plot(y_times,y_amps,'k',linewidth=0.75,zorder=2,
                                                        label=os.path.basename(Y))
        plt.plot(d_times,d_amps,'r',zorder=1,label=os.path.basename(D))
        plt.legend()
        plt.title('{} (black) vs {} (red)'.format(folder_1, folder_2))
        # plt.savefig('{}.png'.format(os.path.basename(Y)))
        plt.show()

if __name__ == "__main__":
    using_obspy()
