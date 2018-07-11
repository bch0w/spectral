"""10.7.18 - Publication level figure generation.

Plot waveforms of all three components for a given GEONET permanent station,
or a temporary RDF station, with preprocessing and filtering set in the script.
Also subplots of determining duration criteria, which is captured using an
amplitude threshold criteria and summing up the time sections (dt) where the
waveform crosses this amplitude threshold.

"""
import os
import sys
import glob
import obspy
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from obspy.core.stream import Stream
from obspy import read, read_inventory, UTCDateTime
from obspy.clients.fdsn import Client

sys.path.append('../modules')
import getdata
import mapmod
from getdata import pathnames
from plotmod import pretty_grids, build_color_dictionary

# ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# global plot parameters
mpl.rcParams['font.size'] = 10
mpl.rcParams['lines.linewidth'] = 1.75

# ============================= RANDOM ONE OFF =================================
def create_obp_nparray():
    """convert the coordinates of obp array into a numpy array for easier hands
    """
    OBPDATA = [['TU13-1', -38, 54.705, 178, 58.573, 3495],
                ['TU13-2', -38, 51.689, 178, 52.993, 2380],
                ['TU13-3', -38, 52.196, 178, 42.740, 1088],
                ['TU13-4', -38, 42.445, 178, 39.617, 1041],
                ['SBPR-1', -38, 43.2603, 178, 53.5904, 2453],
                ['SBPR-2', -38, 50.8454, 178, 52.5146, 2116],
                ['SBPR-3', -38, 53.5952, 178, 45.3122, 1360],
                ['SBPR-4', -38, 54.4714, 178, 58.9590, 3466],
                ['EBPR-1', -38, 44.770, 178, 40.828, 988.6],
                ['EBPR-2', -38, 43.780, 178, 37.134, 1013],
                ['EBPR-3', -38, 41.649, 178, 39.043, 1031],
                ['KU15-1', -38, 54.473, 178, 58.954, 3482.8],
                ['KU15-2', -38, 50.870, 178, 52.405, 2146.6],
                ['KU15-3', -38, 53.260, 178, 45.399, 1363.0],
                ['KU15-4', -38, 42.457, 178, 39.977, 1051.6],
                ['KU15-5', -38, 43.329, 178, 53.970, 2468.6],
                ['KU16-1', -38, 54.437, 178, 59.005, 3473.3],
                ['KU16-2', -38, 50.795, 178, 52.338, 2137.7],
                ['KU16-3', -38, 53.492, 178, 45.327, 1378.5],
                ['KU16-4', -38, 42.676, 178, 39.658, 1047.0],
                ['KU16-5', -38, 43.246, 178, 53.697, 2468.0]]
    
    names,lats,lons,depths = [],[],[],[]
    for obp in OBPDATA:
        names.append(obp[0])
        lat = obp[1] + obp[2]/60
        lon = obp[3] + obp[4]/60
        lats.append(lat)
        lons.append(lon)
        depths.append(obp[4])
        
    dictout = {"NAME":names,"LAT":lats,"LON":lons,"DEPTH":depths}
    np.savez('ITO_OBP_coords.npz',**dictout)
    
def create_lobs_nparray():
    from obspy import read_inventory
    names,lats,lons,depths = [],[],[],[]
    inv = read_inventory('plotYH_HOBITSS_stations_only.xml')
    for net in inv:
        for sta in net:
            names.append(sta.code)
            lats.append(sta.latitude)
            lons.append(sta.longitude)
            depths.append(sta.elevation*-1)
    
    dictout = {"NAME":names,"LAT":lats,"LON":lons,"DEPTH":depths}
    pathout = pathnames()['data'] + 'STATIONXML'
    fullout = os.path.join(pathout,'EBS_LOBS_coords.npz')
    
    np.savez(fullout,**dictout)
            

# ======================== PLOTTING HELPER FUNCTIONS ===========================
def setup_plot():
    """set up plot for waveforms
    """
    f,(ax1,ax2) = plt.subplots(2,sharex=True,sharey=False,figsize=(9,5))
    for ax in [ax1,ax2]:
        pretty_grids(ax)
    f.subplots_adjust(hspace=0.1)
    ax1.set_ylabel('vertical Velocity ($\mu$m/s)')
    ax2.set_ylabel('ground motion ($\mu$m/s)')
    ax2.set_xlabel('time (s)')
    
    return f,ax1,ax2

def plot_polygon(fig):
    # plot low velocity overlay with a polygon
    nw_lat,nw_lon = -39.8764, 178.5385
    ne_lat,ne_lon = -39.0921, 177.1270
    se_lat,se_lon = -37.5630, 178.5248
    sw_lat,sw_lon = -38.3320, 179.9157
    x1,y1 = fig.bmap(nw_lon,nw_lat)
    x2,y2 = fig.bmap(ne_lon,ne_lat)
    x3,y3 = fig.bmap(se_lon,se_lat)
    x4,y4 = fig.bmap(sw_lon,sw_lat)
    poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],
                    facecolor='red',
                    edgecolor='k',
                    linewidth=3,
                    alpha=0.1)
    plt.gca().add_patch(poly)
    
            
def event_beachball_durationmap(eventid,m,anno=False):
    """plt event beachball on figure object, stolen from plotmod.py
    """
    from obspy.imaging.beachball import beach
    
    MT = getdata.get_moment_tensor(eventid)
    eventx,eventy = m(MT['Longitude'],MT['Latitude'])
    FM = [MT['strike1'],MT['dip1'],MT['rake1']]

    b = beach(FM,xy=(eventx,eventy),width=3.5E4,linewidth=1,facecolor='r')
    b.set_zorder(10)
    ax = plt.gca()
    ax.add_collection(b)
    if anno:
        plt.annotate("{}".format(eventid),
                        xy=(eventx,eventy),
                        xytext=(eventx,eventy),
                        fontsize=7,
                        zorder=200,
                        weight='bold')
    
    
def generate_duration_map(event_id):
    """full function to initiate and populate basemap
    """
    map_corner_dict = {"NZ":[-50,-32.5,165,180],
                     "NORTHISLAND":[-42,-34,172,180],
                     "SOUTHISLAND":[-47.5,-40,165,175]}

    f,m = mapmod.initiate_basemap(map_corners=map_corner_dict["NORTHISLAND"],
                                                draw_lines=False)

    # geonet stations
    npzfolder = pathnames()['data'] + 'durationPaper'
    stationfile = '{}_durations.npz'.format(event_id)
    sta = np.load(os.path.join(npzfolder,stationfile))
    names,lats,lons,durs = sta['NAME'],sta['LAT'],sta['LON'],sta['DURATION']
    
    norm = mpl.colors.Normalize(vmin=np.nanmin(durs),vmax=np.nanmax(durs))
    cmap = cm.jet
    colormap = cm.ScalarMappable(norm=norm,cmap=cmap)
    for i,D in enumerate(durs[:75]):
        if np.isnan(D):
            color = 'w'
        else:
            color = colormap.to_rgba(durs[i])
        X,Y = m(lons[i],lats[i])
        m.scatter(X,Y,marker='v',
                        color=color,
                        s=75,
                        edgecolor='k',
                        zorder=5)
        
    
    # event information
    event_beachball_durationmap(event_id,m)
    
    # colorbar
    colormap.set_array(durs)
    cbar = f.colorbar(colormap)
    cbar.set_label('duration (s)',rotation=270,labelpad=8,fontsize=14)
    
    plt.show()
    
    
# ==================================== FUNC ==================================
def get_ito_data(event_id,station="EBPR-1"):
    """grab OBP data from ito
    """
    datapath = pathnames()["mseeds"] + 'ITO'
    fid = "{e}_{s}_vel.mseed".format(e=event_id,s=station)
    filepath = os.path.join(datapath,fid)
    if not os.path.exists:
        print("station event pair doesn't exist")
        sys.exit()
    st = read(filepath)
    
    return st
    
def get_lobs_data(event_id,code):
    """might be able to do this with event_stream?
    """
    excode = "YH.LOBS1"
    st,inv = getdata.fdsn_download(code,event_id,response=True)
    
    return st, inv
    
def get_geonet_data(event_id,code):
    # GATHER
    st,inv,cat = getdata.event_stream(code=code,event_id=event_id,startpad=15)
    event = cat[0]
    
    return st, inv, cat
    
def preprocess(st_in,inv=None):
    """grab data using external functions
    """
    st = st_in.copy()
    
    # PREPROCESS and CONVERT to cm/s
    st.detrend('linear')
    st.detrend('demean')
    st.taper(max_percentage=0.05)
    if inv:
        st.attach_response(inventories=inv)
        st.remove_response(output='VEL',water_level=60)
    
    st.detrend('linear')
    st.detrend('demean')
    st.taper(max_percentage=0.05)

    for tr in st:
        tr.data *= 1E6

    # FILTER
    tmin,tmax = 10,100
    freqmin,freqmax = 1/tmax,1/tmin
    st.filter("bandpass",freqmin=freqmin,freqmax=freqmax,corners=3)
    st.detrend('linear')
    st.detrend('demean')
    st.taper(max_percentage=0.05)

    return st

def amplitude_threshold(st,threshold_percentage):
    """determine amplitudes over a given threshold and count the time length
    for amount of trace over this threshold. return arrays for plotting
    """
    # for geonet data
    if len(st) == 3:
        groundmotion = np.sqrt(st[0].data**2 + st[1].data**2 + st[2].data**2)
        vertical = st.select(component='Z')[0].data
        
    # for obp data with only 1 comp
    elif len(st) == 1:
        groundmotion = np.sqrt(st[0].data**2)
        vertical = st[0].data
    
    peakamplitude = groundmotion.max()
    threshold = peakamplitude * threshold_percentage

    # array of values of threshold
    whereover = np.where(groundmotion >= threshold)
    
    # groundmotion over threshold for plotting
    groundmotion_over = np.copy(groundmotion)
    groundmotion_over[groundmotion_over<threshold] = np.nan
        
    return vertical,groundmotion,whereover,groundmotion_over,threshold
    
    
# ================================== PROCESSING  ===============================
def process_and_plot_waveforms(event_id,code,threshold_choice=0.2,choice="GN",
                                                        show=False,save=False):
    """create first section of composite, two waveforms
    """
    f,ax1,ax2 = setup_plot()
    
    # process 
    if choice == "GN":
        st_raw,inv,event = get_geonet_data(event_id,code)
        st = preprocess(st_raw,inv)
        code = code.split('.')[1]
    elif choice == "ITO":
        st_raw = get_ito_data(event_id,code)
        _,_,event = get_geonet_data(event_id,'NZ.PUZ.*.HH?')
        st = preprocess(st_raw)
        st.resample(20)
        for ax in [ax1,ax2]:
            ax.set_ylabel('velocity ($\mu$m/s)')
    elif choice == "YH":
        st_raw,inv  = get_lobs_data(event_id,code)
        _,_,event = get_geonet_data(event_id,'NZ.PUZ.*.HH?')
        st = preprocess(st_raw,inv)
        code = code.split('.')[1]

    vertical,groundmotion,whereover,groundmotion_over,threshold = \
                                        amplitude_threshold(st,threshold_choice)

    # time axes
    stats= st[0].stats
    t = np.linspace(0,stats.endtime-stats.starttime,stats.npts)
    duration = len(whereover[0]) / stats.sampling_rate

    # plot
    ax1.set_title(code)
    ax1.plot(t,vertical,'k')
    ax2.plot(t,groundmotion,'k',zorder=3)
    ax2.plot(t,groundmotion_over,'r',zorder=4,
                label='Duration = {}s'.format(int(duration)))
    ax2.axhline(xmin=t[0],xmax=t[-1],y=threshold,zorder=2,
                color='gray',linestyle='-.',linewidth=1.5,
                label='20% peak amplitude') 
    ax2.legend()   
    # set bounds
    for ax in [ax1,ax2]:
        ax.set_xlim([t[0],t[-1]])
    
    if save:
        outputfolder = pathnames()['spectralplots'] + 'durations/PAPEROUT'
        fid = 'wav_{e}_{c}.png'.format(e=event_id,c=code)
        fidout = os.path.join(outputfolder,fid)
        plt.savefig(fidout,dpi=100)
    if show:
        plt.show()

    plt.close()
    return duration

# ================================== RUN SCRIPTS ==============================
def run_waveform_plotter():
    """process_and_plot_waveforms
    """
    event_ids = ['2014p864702','2015p822263','2014p240655']
    
    # GEONET DATA
    code = 'NZ.PUZ.*.HH?'
    duration = process_and_plot_waveforms(event_ids[0],code,choice='GN')
    
    # ITO-SAN DATA
    itocode = 'EBPR-2'
    duration = process_and_plot_waveforms(event_ids[0],itocode,choice='ITO')
    
def loop_waveform_plotter_geonet(event_id):
    """generate output file to be read in for mapping
    """
    # collect station information
    npzfolder = pathnames()['data'] + 'STATIONXML'
    stationfile = 'NZ_BB_coords.npz'
    sta = np.load(os.path.join(npzfolder,stationfile))
    names,lats,lons = sta['NAME'],sta['LAT'],sta['LON']
    
    durations = []
    for sta in names:
        try:
            code = 'NZ.{}.*.HH?'.format(sta)
            duration = process_and_plot_waveforms(event_id,code,choice='GN',
                                                        show=False,save=False)
            durations.append(duration)
        except Exception as e:
            print(sta)
            durations.append(np.nan)
            plt.close()
            continue
    
    # save into np file for quick access
    dictout = {"NAME":names,"DURATION":durations,"LAT":lats,"LON":lons}
    pathout = pathnames()['data'] + 'durationPaper'
    fidout = "{}_durations".format(event_id)
    allout = os.path.join(pathout,fidout)
    np.savez(allout,**dictout)
        
        
def loop_waveform_plotter_itodata(event_id):
    """generate output file to be read in for mapping, to be run after geonet 
    data to access an already created npz file
    """    
    # collect event file information
    mseedpath = pathnames()['mseeds'] + 'ITO/{}*'.format(event_id)
    eventfiles = glob.glob(mseedpath)
    stationnames = []
    for fid in eventfiles:
        stationnames.append(fid.split('_')[1])
    
    # collect station information
    coordsfile = pathnames()['data'] +'STATIONXML/ITO_OBP_coords.npz'
    ito = np.load(coordsfile)
    itonames,itolats,itolons = ito['NAME'],ito['LAT'],ito['LON']
    
    newnames,newlats,newlons,durations = [],[],[],[]
    for code in stationnames:
        try:
            itoind = np.where(itonames==code)[0][0]
            newnames.append(itonames[itoind])
            newlats.append(itolats[itoind])
            newlons.append(itolons[itoind])
            duration = process_and_plot_waveforms(event_id,code,choice='ITO',
                                                        show=False,save=False)
            durations.append(duration)
        except Exception as e:
            print(sta)
            durations.append(np.nan)
            plt.close()
            continue

    # collect station information
    pathout = pathnames()['data'] + 'durationPaper'
    stationfile = '{}_durations.npz'.format(event_id)
    allout = os.path.join(pathout,stationfile)
    sta = np.load(allout)
    
    names,lats,lons,durs = sta['NAME'],sta['LAT'],sta['LON'],sta['DURATION']
    names = np.concatenate((names,np.array(newnames)))
    lats = np.concatenate((lats,np.array(newlats)))
    lons = np.concatenate((lons,np.array(newlons)))
    durs = np.concatenate((durs,np.array(durations)))
    
    dictout = {"NAME":names,"DURATION":durs,"LAT":lats,"LON":lons}
    np.savez(allout,**dictout)

def loop_waveform_plotter_lobsdata(event_id):
    """generate output file to be read in for mapping, to be run after geonet 
    data to access an already created npz file
    """    
    # collect station information
    coordsfile = pathnames()['data'] +'STATIONXML/EBS_LOBS_coords.npz'
    lobs = np.load(coordsfile)
    lobsnames,lobslats,lobslons = lobs['NAME'],lobs['LAT'],lobs['LON']
    
    newnames,newlats,newlons,durations = [],[],[],[]
    for sta in lobsnames:
        try:
            code = "YH.{}.*.?H?".format(sta)
            lobsind = np.where(lobsnames==sta)[0][0]
            newnames.append(lobsnames[lobsind])
            newlats.append(lobslats[lobsind])
            newlons.append(lobslons[lobsind])
            duration = process_and_plot_waveforms(event_id,code,choice='YH',
                                                        show=False,save=True)
            durations.append(duration)
        except Exception as e:
            print(sta)
            durations.append(np.nan)
            plt.close()
            continue
    
    import ipdb;ipdb.set_trace()
    # collect station information
    pathout = pathnames()['data'] + 'durationPaper'
    stationfile = '{}_durations.npz'.format(event_id)
    allout = os.path.join(pathout,stationfile)
    sta = np.load(allout)
    
    names,lats,lons,durs = sta['NAME'],sta['LAT'],sta['LON'],sta['DURATION']
    names = np.concatenate((names,np.array(newnames)))
    lats = np.concatenate((lats,np.array(newlats)))
    lons = np.concatenate((lons,np.array(newlons)))
    durs = np.concatenate((durs,np.array(durations)))
    
    dictout = {"NAME":names,"DURATION":durs,"LAT":lats,"LON":lons}
    np.savez(allout,**dictout)
        
if __name__ == "__main__":
    event_id = '2014p864702'
    # loop_waveform_plotter_geonet(event_id)
    # loop_waveform_plotter_itodata(event_id)
    # loop_waveform_plotter_lobsdata(event_id)
    generate_duration_map(event_id)

