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
import matplotlib.pyplot as plt

mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1

# ================================= HELPERS ====================================
def setup_plot(number_of_plots):
    """dynamically set up plots according to number of files, modified from
    magnifytremors.py script
    """
    f = plt.figure(figsize=(11.69,8.27),dpi=100)
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
    

# =================================== FUNC ====================================
def waveform_getter(code,mesh,bounds=None):
    """grab mseed data from synthetics folder, filter if necessary
    available inputs for mesh: ETOPO1, SRTM15P, SRTM30P
    """
    # initialization
    net,sta,loc,cha = code.split('.')
    datapath = join(pathnames()['syns'],'2018130600_meshTest',mesh,'')
    datafiles = glob.glob(join(datapath,'{n}.{s}*.mseed'.format(n=net,s=sta)))
    if not datafiles:
        print('No data files for {n}.{s} in mesh {m} found'.format(n=net,
                                                                   s=sta,
                                                                   m=mesh))
        return None
    
    # read in datafiles
    st = Stream()
    for d in datafiles:
        st += read(d)
    
    # filter if necessary
    if bounds:
        st.filter('bandpass',freqmin=bounds[0],freqmax=bounds[1])
    
    return st
    
def waveform_plotter(code,bounds,show=True,save=False):
    """call waveform_getter for the same station in each mesh and plot the 
    waveforms superimposed to see the differences
    """
    net,sta,loc,cha = code.split('.')

    f,axes = setup_plot(3)
    for c,mesh in zip(['r','g','b'],['ETOPO1','SRTM15P','SRTM30P']):
        st = waveform_getter(code,mesh,bounds)
        for ax,comp in zip(axes,['N','E','Z']):
            st_comp = st.select(component=comp)
            plotlabel = "{m} {n}.{s}".format(m=mesh,n=net,s=sta)
            ax.plot(st_comp.data,c=c,label=plotlabel)
            
    for ax,comp in zip(axes,['N','E','Z']):
        ax.set_ylabel('{} disp [m]'.format(comp))
        ax.legend(prop={"size":5})
        
    axes[2].set_xlabel('Time [s]')
    plt.title('Mesh Waveform Comparison')
    
    if show:
        plt.show()
    if save:
        figtitle = "meshcompare_{n}.{s}.png".format(n=net,s=sta)
        figfolder = join(pathnames()['kupeplots'],"meshcompare",figtitle)
        plt.savefig(figfolder,dpi=100)

def mesh_mapper_2d(mesh,show=True,save=False):
    """plot an xyz topography file as a 3d figure. lags out hard because our 
    topography data is so large and high resolution
    """    
    meshes = {"ETOPO1":"topo_ETOPO1_utm60H_ismooth0_553_622_1000m_surf.xyz",
              "SRTM30P":"topo_SRTM30P_utm60H_ismooth0_553_622_1000m_surf.xyz",
              "SRTM15P":"topo_SRTM15P_utm60H_ismooth0_553_622_1000m_surf.xyz"
              }
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
    plt.colorbar()
    
    if show:
        plt.show()
    if save:
        figtitle = "mesh_{m}.png".format(m=mesh)
        figfolder = join(pathnames()['kupeplots'],"meshcompare",figtitle)
        plt.savefig(figfolder,dpi=400) 
    
    return f
       
def event_station_plotter(event,station):
    """given an event and a station, plot them on a map
    """
    # find station coordinates
    nz_bb_path = pathnames()['DATA'] + 'STATIONXML/nz_BB_coords.npz'
    nz_bb = np.load(nz_bb_path)
    nz_bb_names = nz_bb['names']
    
    bb_ind = np.where(nz_bb_names==station)
    sta_lat = nz_bb['lats'][bb_ind][0][0]
    sta_lon = nz_bb['lons'][bb_ind][0][0]
    
    # find event coordinates
    
        
def mesh_differencer_2d(mesh1,mesh2,show=True,save=False):
    """difference the elevation of mesh1 and mesh2
    """
    meshes = {"ETOPO1":"topo_ETOPO1_utm60H_ismooth0_553_622_1000m_surf.xyz",
              "SRTM30P":"topo_SRTM30P_utm60H_ismooth0_553_622_1000m_surf.xyz",
              "SRTM15P":"topo_SRTM15P_utm60H_ismooth0_553_622_1000m_surf.xyz"
              }
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

if __name__ == "__main__":
    mesh_differencer_2d('ETOPO1','SRTM15P',show=True,save=False)
    # mesh_mapper_2d('SRTM30P',show=False,save=True)
        
    
    
    
    
    
    
    
    
    
    
    