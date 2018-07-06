"""in generating meshes with different topography files at different resolutions
we need a way to quantify the different meshes and determine what is different
in the resulting waveforms from each mesh, how the topography files affect
the surface waves, and how the mesh edges affect stations that are very close to
the boundaries. this script should hopefully encompass all the tools to make
that possible
++worklist
-map plotter with topography and bounds as input
-mesh differencer, to produce maps that show where meshes vary
-combined plotter that will put all these figures together
"""
import os
import sys
import glob
import math
sys.path.append('../modules/')
import numpy as np
from obspy import read, Stream
from os.path import join

from getdata import pathnames
from plotmod import pretty_grids

import matplotlib as mpl
# mpl.use('TkAgg')
from matplotlib import gridspec
import matplotlib.pyplot as plt

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.25

# ================================= HELPERS ====================================
def setup_plot(number_of_plots):
    """dynamically set up plots according to number of files, modified from
    magnifytremors.py script
    """
    f = plt.figure(figsize=(11.69,8.27),dpi=75)
    nrows,ncols=number_of_plots,1
    height_ratios = [1] * (number_of_plots)

    gs = gridspec.GridSpec(nrows,ncols,height_ratios=height_ratios,hspace=0)

    # create gridspec subplots, sharex with the first axis
    axes,twaxes = [],[]
    for i in range(number_of_plots):
        if i == 0:
            ax = plt.subplot(gs[i])
        else:
            ax = plt.subplot(gs[i],sharex=axes[0])

        pretty_grids(ax)
        axes.append(ax)

    # remove x tick labels except for last axis
    for ax in axes[0:-1]:
        plt.setp(ax.get_xticklabels(),visible=False)

    return f,axes

def lonlat_utm(lon_or_x,lat_or_y,zone=60,inverse=False):
    """!!!taken from meshGenTools
    convert latitude and longitude coordinates to UTM projection

    :type lon_or_x: float or int
    :param lon_or_x: longitude value in WGS84 or X in UTM-'zone' projection
    :type lat_or_y: float or int
    :param lat_or_y: latude value in WGS84 or Y in UTM-'zone' projection
    :type zone: int
    :param zone: UTM zone for conversion from WGS84
    :type inverse: bool
    :param inverse: if inverse == False, latlon => UTM, vice versa.
    """
    from pyproj import Proj
    projstr = ("+proj=utm +zone=60, +south +ellps=WGS84 +datum=WGS84"
                                                        " +units=m +no_defs")
    myProj = Proj(projstr)
    x_or_lon,y_or_lat = myProj(lon_or_x,lat_or_y,inverse=inverse)

    return x_or_lon, y_or_lat


# =================================== MESH ====================================
def mesh_differencer_2d(mesh1,mesh2,show=True,save=False):
    """difference the elevation of mesh1 and mesh2
    """

    # MESH 1
    filepath1 = pathnames()['xyz'] + meshes[mesh1]
    # filepath1=('topo_SRTM30P_utm60H_ismooth0_553_622_1000m_surf.xyz')
    data1 = np.genfromtxt(filepath1)
    x1 = data1[:,0]
    y1 = data1[:,1]
    z1 = data1[:,2]

    # MESH 2
    filepath2 = pathnames()['xyz'] + meshes[mesh2]
    # filepath2=('topo_SRTM15P_utm60H_ismooth0_553_622_1000m_surf.xyz')
    data2 = np.genfromtxt(filepath2)
    x2 = data2[:,0]
    y2 = data2[:,1]
    z2 = data2[:,2]

    zdiff = z1-z2

    # plot it
    f,ax = plt.subplots()
    plt.scatter(x1,y1,c=zdiff,cmap='seismic')

    # formatting
    plt.ticklabel_format(style='sci',
                         axis='both',
                         scilimits=(0,0))
    plt.xlim([min(x1),max(x1)])
    plt.ylim([min(y1),max(y1)])
    plt.title("Mesh difference between {m1} and {m2}".format(m1=mesh1,
                                                             m2=mesh2)
                                                             )
    globalmax = math.ceil(max([min(zdiff),max(zdiff)]))
    norm = mpl.colors.Normalize(vmin=min(zdiff),vmax=max(zdiff))
    cbar = plt.colorbar()#boundaries=
                                # np.linspace(min(zdiff),max(zdiff),201))
    cbar.ax.set_ylabel('elevation difference [m]', rotation=90)
    if show:
        plt.show()
    if save:
        figtitle = "meshdiff_{m1}_{m2}.png".format(m1=mesh1,m2=mesh2)
        figfolder = join(pathnames()['kupeplots'],"meshcompare",figtitle)
        plt.savefig(figfolder,dpi=400)


def mesh_mapper_3d():
    """plot an xyz topography file as a 3d figure. lags out hard because our
    topography data is so large and high resolution
    """
    from matplotlib import cm
    from matplotlib.mlab import griddata
    from mpl_toolkits.mplot3d import Axes3D

    filepath = (pathnames()['xyz'] +
                        'topo_SRTM30P_utm60H_ismooth0_553_622_1000m_surf.xyz')
    data = np.genfromtxt(filepath)

    x = data[:343966,0]
    y = data[:343966,1]
    z = data[:343966,2]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    xi = np.linspace(min(x), max(x))
    yi = np.linspace(min(y), max(y))

    X, Y = np.meshgrid(xi, yi)
    Z = griddata(x, y, z, xi, yi)

    surf = ax.plot_surface(X, Y, Z, rstride=5, cstride=5, cmap=cm.jet,
                           linewidth=1)

    ax.set_zlim3d(np.min(Z), np.max(Z))
    fig.colorbar(surf)

    plt.show()

def mesh_mapper_2d(mesh,show=True,save=False):
    """plot an xyz topography file as a 3d figure. lags out hard because our
    topography data is so large and high resolution
    """
    filepath = pathnames()['xyz'] + meshes[mesh]

    data = np.genfromtxt(filepath)

    x = data[:,0]
    y = data[:,1]
    z = data[:,2]

    f = plt.scatter(x,y,c=z,cmap='inferno')
    plt.ticklabel_format(style='sci',
                         axis='both',
                         scilimits=(0,0))
    plt.xlim([min(x),max(x)])
    plt.ylim([min(y),max(y)])
    plt.title(mesh)
    cbar = plt.colorbar(format='%.0e')
    cbar.ax.set_ylabel('elevation [m]', rotation=90)

    if save:
        figtitle = "mesh_{m}.png".format(m=mesh)
        figfolder = join(pathnames()['kupeplots'],"meshcompare",figtitle)
        plt.savefig(figfolder,dpi=400)
    if show:
        plt.show()

    return f

# ================================ WAVEFORM ====================================
def waveform_getter(code,mesh,bounds=None):
    """grab mseed data from synthetics folder, filter if necessary
    available inputs for mesh: ETOPO1, SRTM15P, SRTM30P
    """
    # initialization
    net,sta,loc,cha = code.split('.')
    datapath = join(pathnames()['syns'],'2018p130600_meshTest',mesh,'')
    datafiles = glob.glob(join(datapath,'{n}.{s}*.mseed'.format(n=net,s=sta)))

    # data checkstop
    if not datafiles:
        print('No data files for {n}.{s} in mesh {m} found'.format(n=net,
                                                                   s=sta,
                                                                   m=mesh))
        return None

    # read in datafiles
    st = Stream()
    for d in datafiles:
        st += read(d)

    # filter if necessary, convert bounds from seconds to frequency
    if bounds:
        freqmin = 1/bounds[1]
        freqmax = 1/bounds[0]
        st.filter('bandpass',freqmin=freqmin,freqmax=freqmax)

    return st

def waveform_plotter(code,bounds,meshlist,show=True,save=False):
    """call waveform_getter for the same station in each mesh and plot the
    waveforms superimposed to see the differences
    """
    net,sta,loc,cha = code.split('.')

    colordict = {"ETOPO1_553_622":'r',
                 "ETOPO1_496_566":'c',
                 "SRTM30P_552_620":'g',
                 "SRTM15P":'b'
                 }

    # dynamically plot waveforms based on meshlist
    f,axes = setup_plot(3)
    for i,mesh in enumerate(meshlist):
        st = waveform_getter(code,mesh,bounds)
        if i == 0:
            stats = st[0].stats
            t = np.linspace(0,stats.endtime-stats.starttime,stats.npts)
        for ax,comp in zip(axes,['N','E','Z']):
            st_comp = st.select(component=comp)
            plotlabel = "{m} {n}.{s}".format(m=mesh,n=net,s=sta)
            ax.plot(t,st_comp[0].data,c=colordict[mesh],label=plotlabel)

    # post-formatting
    for ax,comp in zip(axes,['N','E','Z']):
        ax.set_ylabel('{} disp [m]'.format(comp))
        ax.legend(prop={"size":10})
    axes[-1].set_xlabel('Time [s]')
    plt.xlim([0,300])
    plt.sca(axes[0])
    plt.title('Mesh Waveform Comparison [{b0}-{b1}s]'.format(b0=bounds[0],
                                                           b1=bounds[1]))

    if save:
        figtitle = "meshwaveforms_{n}.{s}_{b0}-{b1}.png".format(n=net,s=sta,
                                                                b0=bounds[0],
                                                                b1=bounds[1])
        figfolder = join(pathnames()['kupeplots'],"meshcompare",figtitle)
        plt.savefig(figfolder,dpi=100)

    if show:
        plt.show()

    return f,axes


def event_station_plotter(station,mesh='SRTM30P',
                                    event="2018p130600",save=False,show=True):
    """given an event and a station, plot them on a map
    """
    # find station coordinates
    nz_bb_path = pathnames()['data'] + 'STATIONXML/nz_BB_coords.npz'
    nz_bb = np.load(nz_bb_path)
    nz_bb_names = nz_bb['NAME']

    bb_ind = np.where(nz_bb_names==station)[0][0]
    sta_lat = nz_bb['LAT'][bb_ind]
    sta_lon = nz_bb['LON'][bb_ind]

    # static event coordinates - but could get from CMTSOLUTIONS easily
    ev_lat = -39.9490
    ev_lon = 176.2995

    # convert coordinates to UTM60
    sta_x,sta_y = lonlat_utm(sta_lon,sta_lat)
    ev_x,ev_y = lonlat_utm(ev_lon,ev_lat)

    # call mesh plotter and superimpose station and event
    f = mesh_mapper_2d(mesh,show=False,save=False)
    # plt.figure(f)
    plt.plot([sta_x,ev_x],[sta_y,ev_y],'k',zorder=4)
    plt.scatter(sta_x,sta_y,s=22.5,c='b',marker='v',
                                        edgecolors='k',linewidths=1,zorder=5)
    plt.annotate(s="NZ.{}".format(station),xy=(sta_x,sta_y),zorder=6,color='w')
    plt.scatter(ev_x,ev_y,s=22.5,c='r',marker='o',
                                        edgecolors='k',linewidths=1,zorder=5)
    plt.annotate(s='2018p130600',xy=(ev_x,ev_y),zorder=6,color='w')

    plt.xlabel('X (UTM-60)')
    plt.ylabel('Y (UTM-60)')

    if save:
        figtitle = "evsta_{s}.png".format(s=station)
        figfolder = join(pathnames()['kupeplots'],"meshcompare",figtitle)
        plt.savefig(figfolder,dpi=100)
    if show:
        plt.show()

    return f



# ============================== FUNC CALLS ====================================
def create_composite():
    """waveform comparisons of synthetics with a map showing source-receiver
    pair and connecting lines, potential option to plot onto real mesh (might
    be intensive)
    """

    # meshlist = ["ETOPO1_496_566","SRTM30P_552_620","ETOPO1_553_622",]
    meshlist = ["SRTM30P_552_620","ETOPO1_553_622",]

    boundlist=[[3,30],[10,30]]

    # create station list
    synth_list = (pathnames()['syns'] +
                            '2018p130600_meshTest/ETOPO1_496_566/*.mseed')
    stations = glob.glob(synth_list)
    station_list = []
    for s in stations:
        _,S,_,_,_ = s.split('.')
        station_list.append(S)

    station_list = list(set(station_list))
    station_list = ["PUZ","MWZ","NNZ","TOZ"]

    for sta in station_list:
        # try:
        code = "NZ.{}..HHZ".format(sta)
        f_es = event_station_plotter(sta,show=False,save=True)
        plt.close()
        for bounds in boundlist:
            f_wf,axes = waveform_plotter(code,bounds,

                                            meshlist,show=False,save=True)
            plt.close()
        # except Exception as e:
        #     print(e)
        #     continue




if __name__ == "__main__":
    # multi-use mesh dictionary
    global meshes
    meshes = {
    "ETOPO1_553_622":"topo_ETOPO1_utm60H_ismooth0_553_622_1000m_surf.xyz",
    "SRTM30P":"topo_SRTM30P_utm60H_ismooth0_553_622_1000m_surf.xyz",
    "SRTM30P_552":"topo_SRTM30P_utm60H_ismooth0_552_620surf.xyz",
    "SRTM15P":"topo_SRTM15P_utm60H_ismooth0_553_622_1000m_surf.xyz",
    "ETOPO1_496_566":"topo_NZ_BC_utm60H_ismooth0_496_566surf.xyz"
              }

    create_composite()
