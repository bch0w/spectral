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
from procmod import myround
from plotmod import pretty_grids, build_color_dictionary

# ignore warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# global plot parameters
mpl.rcParams['font.size'] = 11
mpl.rcParams['lines.linewidth'] = 1.75

# ============================= RANDOM ONE OFF =================================
def convert_ITODATA_to_mseed():
    """ito-san gave us sac data with no timestamps and in acceleration.
    set proper time stamp, convert to velocity and save as mseed
    """
    import glob
    from obspy import read, UTCDateTime


    timedict = {'2014p864702': UTCDateTime('2014-11-16T22:33:17'),
                '2014p051675': UTCDateTime('2014-01-20T02:52:45'),
                '2015p768477': UTCDateTime('2015-10-12T08:05:01'),
                '2014p715167': UTCDateTime('2014-09-22T14:41:22'),
                '2014p240655': UTCDateTime('2014-03-31T01:01:19'),
                '2016p859524': UTCDateTime('2016-11-14T00:34:22')
                }

    sacfiles = glob.glob('*.s')
    for sf in sacfiles:
        st = read(sf)
        event = sf.split('_')[2].split('.')[0]
        st[0].stats.starttime = timedict[event]
        newname = "{}_{}_vel.mseed".format(event,sf.split('_')[0])
        st.differentiate()
        for tr in st:
            tr.stats.station = sf.split('_')[0]
            tr.stats.network = "YH"
        st.write(newname,format="MSEED")

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
    st,inv = getdata.fdsn_download(code,event_id,response=True)

    return st, inv

def get_geonet_data(event_id,code):
    # GATHER
    st,inv,cat = getdata.event_stream(code=code,event_id=event_id,startpad=15)
    event = cat[0]

    return st, inv, event

def preprocess(st_in,inv=None,event=None):
    """grab data using external functions
    """
    st = st_in.copy()

    # trim to event origin time
    if event:
        origintime = event.origins[0].time
        st.trim(origintime,origintime+1000)
        
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

    # FILTER, global bounds
    freqmin,freqmax = 1/bounds[1],1/bounds[0]
    st.filter("bandpass",freqmin=freqmin,freqmax=freqmax,corners=3)
    st.detrend('linear')
    st.detrend('demean')
    st.taper(max_percentage=0.05)
    

    return st

def amplitude_threshold(st,threshold_percentage):
    """determine amplitudes over a given threshold and count the time length
    for amount of trace over this threshold. return arrays for plotting
    """
    # for 3 component data
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

# ======================== PLOTTING HELPER FUNCTIONS ===========================
def setup_plot():
    """set up plot for waveforms
    """
    f,(ax1,ax2) = plt.subplots(2,sharex=True,sharey=False,figsize=(9,5))
    for ax in [ax1,ax2]:
        pretty_grids(ax)
    f.subplots_adjust(hspace=0.1)
    ax1.set_ylabel('vertical velocity ($\mu$m/s)')
    ax2.set_ylabel('ground motion ($\mu$m/s)')
    ax2.set_xlabel('time (s)')

    return f,ax1,ax2
    
def trace_trench(m):
    """trace the hikurangi trench on a basemap object
    """
    trenchcoordspath = pathnames()['data'] + \
                                        'DURATIONS/hikurangiTrenchCoords.npz'
    trenchcoords = np.load(trenchcoordspath)
    lats = trenchcoords['LAT']
    lons = trenchcoords['LON']
    x,y = m(lons,lats)
    
    # interpolate points to make a smoother curve
    xprime = np.flip(x,axis=0)
    yprime = np.flip(y,axis=0)
    xprimenew = np.linspace(x.min(),x.max(),100)
    yprimenew = np.interp(xprimenew,xprime,yprime)

    m.plot(xprimenew,yprimenew,'--',linewidth=0.9,color='k',zorder=2) 

def event_beachball_durationmap(eventid,m,anno=False):
    """plt event beachball on figure object, stolen from plotmod.py
    """
    from obspy.imaging.beachball import beach

    MT = getdata.get_moment_tensor(eventid)
    eventx,eventy = m(MT['Longitude'],MT['Latitude'])
    FM = [MT['strike2'],MT['dip2'],MT['rake2']]

    b = beach(FM,xy=(eventx,eventy),width=3.2E4,linewidth=1,facecolor='r')
    b.set_zorder(10)
    ax = plt.gca()
    ax.add_collection(b)
    if anno:
        plt.annotate("M{m}\n{e}".format(m=MT['Mw'],e=eventid),
                     xy=(eventx,eventy),
                     xytext=(eventx,eventy),
                     fontsize=10,
                     zorder=200,
                     weight='bold')


def generate_duration_map(corners,event_id):
    """full function to initiate and populate basemap
    """
    manual_ignore = ["LOBS1","LOBS4",]

    f,m = mapmod.initiate_basemap(map_corners=corners,draw_lines=False)

    # load in duration values
    npzfolder = pathnames()['data'] + 'DURATIONS'
    stationfile = '{e}_{b0}_{b1}_durations.npz'.format(e=event_id,
                                                       b0=bounds[0],
                                                       b1=bounds[1]
                                                       )
    sta = np.load(os.path.join(npzfolder,stationfile))

    # manual ignorance
    names,lats,lons,durs = sta['NAME'],sta['LAT'],sta['LON'],sta['DURATION']
    for MI in manual_ignore:
        ind = np.where(names==MI)[0][0]
        names = np.delete(names,ind)
        lats = np.delete(lats,ind)
        lons = np.delete(lons,ind)
        durs = np.delete(durs,ind)

    # set colormap to the values of duration
    vmax = myround(np.nanmax(durs),base=50,choice='up')
    vmin = myround(np.nanmin(durs),base=50,choice='down')
    norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    cmap = cm.jet
    colormap = cm.ScalarMappable(norm=norm,cmap=cmap)

    # plot individual stations with color and zorder set by duration length
    for i,D in enumerate(durs):
        if np.isnan(D):
            color = 'w'
            D = 10
        else:
            color = colormap.to_rgba(durs[i])
        X,Y = m(lons[i],lats[i])
        m.scatter(X,Y,marker='v',
                        color=color,
                        s=75,
                        edgecolor='k',
                        zorder=int(D))
        if (names[i] == 'OPRZ') or (names[i] == 'LOBS6'):
            plt.annotate(names[i],
                         xy=(X,Y),
                         fontsize=10,
                         zorder=1000,
                         weight='bold')

    # additional map objects
    trace_trench(m)
    event_beachball_durationmap(event_id,m,anno=True)
    m.drawmapscale(179,-41.75, 179,-41.75, 100,yoffset=0.01*(m.ymax-m.ymin))

    # colorbar
    colormap.set_array(durs)
    cbar = f.colorbar(colormap)
    cbar.set_label('duration (s)',rotation=270,labelpad=17,fontsize=11.5)

    plt.title(event_id)
    plt.show()

# ================================== PROCPLOT =================================
def true_if_outside_bounds(lat,lon,corners):
    """if station falls outside the map bounds, don't process. set up for new
    zealand coordinates
    """    
    lat_bot,lat_top,lon_left,lon_right = corners
    
    if (lon < lon_left) or (lon > lon_right):
        return True
    elif (lat < lat_bot) or (lat > lat_top):
        return True
    else:
        return False

def process_and_plot_waveforms(event_id,code,threshold_choice=0.2,choice="GN",
                                                        show=False,save=False):
    """create first section of composite, two waveforms
    """
    f,ax1,ax2 = setup_plot()

    # PROCESS by network
    if choice == "GN":
        st_raw,inv,event = get_geonet_data(event_id,code)
        st = preprocess(st_raw,inv,event)
        code = code.split('.')[1]
        x0,x1 = 0,1000
    elif choice == "ITO":
        st_raw = get_ito_data(event_id,code)
        _,_,event = get_geonet_data(event_id,'NZ.PUZ.*.HH?')
        st = preprocess(st_raw,event=event)
        st.resample(20)
        x0,x1 = 1000,2000
        for ax in [ax1,ax2]:
            ax.set_ylabel('velocity ($\mu$m/s)')
    elif choice == "YH":
        st_raw,inv  = get_lobs_data(event_id,code)
        _,_,event = get_geonet_data(event_id,'NZ.PUZ.*.HH?')
        st = preprocess(st_raw,inv,event)
        code = code.split('.')[1]
        x0,x1 = 0,1000

    vertical,groundmotion,whereover,groundmotion_over,threshold = \
                                        amplitude_threshold(st,threshold_choice)

    # calculate duration from counted samples
    stats= st[0].stats
    duration = len(whereover[0]) / stats.sampling_rate
    if (show==False) and (save==False):
        plt.close()
        return duration
        
    # PLOT WAVEFORMS
    t = np.linspace(0,stats.endtime-stats.starttime,stats.npts)
    ax1.plot(t,vertical,'k')
    ax2.plot(t,groundmotion,'k',zorder=3)
    ax2.plot(t,groundmotion_over,'r',zorder=4,
                label='Duration = {}s'.format(int(duration)))
    ax2.axhline(xmin=t[0],xmax=t[-1],y=threshold,zorder=2,
                color='gray',linestyle='-.',linewidth=1.5,
                label='20% peak amplitude')
    
    # formatting
    ax1.set_title(code + ' [{b0}-{b1}s]'.format(b0=bounds[0],b1=bounds[1]))
    ax2.legend()
    for ax in [ax1,ax2]:
        ax.set_xlim([x0,x1])

    # checksave
    if save:
        outputfolder = pathnames()['spectralplots'] + 'durations/PAPEROUT'
        fid = 'wav_{e}_{c}_{b0}_{b1}.png'.format(e=event_id,c=code,
                                                 b0=bounds[0],b1=bounds[1])
        fidout = os.path.join(outputfolder,fid)
        plt.savefig(fidout,dpi=100)
    if show:
        plt.show()

    plt.close()
    
    return duration

# ================================== RUN SCRIPTS ==============================
def loop_waveform_plotter(event_id,corners,choice='GN',show=True,save=False):
    """generate output file to be
            ax.set_xlim([400,1400]) read in for mapping
    """
    filedict = {"GN":"NZ_BB_coords.npz",
                "ITO":"ITO_OBP_coords.npz",
                "YH":"LOBS_coords.npz"}
    stationfile = filedict[choice]

    # collect station information
    npzfolder = pathnames()['data'] + 'STATIONXML'
    sta = np.load(os.path.join(npzfolder,stationfile))
    names,lats,lons = sta['NAME'],sta['LAT'],sta['LON']

    # numpy array for saving duration information
    pathout = pathnames()['data'] + 'DURATIONS'
    stationfile = '{e}_{b0}_{b1}_durations.npz'.format(e=event_id,
                                                       b0=bounds[0],
                                                       b1=bounds[1]
                                                       )
    allout = os.path.join(pathout,stationfile)

    # geonet data starts off the numpy list
    if choice == "GN":
        durations = []
        for i,sta in enumerate(names):
            try:
                if true_if_outside_bounds(lats[i],lons[i],corners):
                    print(sta,'outside bounds')
                    continue
                
                code = 'NZ.{}.*.HH?'.format(sta)
                duration = process_and_plot_waveforms(event_id,code,
                                                        choice=choice,
                                                        show=show,
                                                        save=save)
                durations.append(duration)
            except Exception as e:
                print('\t'+sta)
                durations.append(np.nan)
                plt.close()
                continue
        dictout = {"NAME":names,"DURATION":durations,"LAT":lats,"LON":lons}


    # ITO and LOBS data will append to the numpy list
    else:
        # import numpy array to be appended to
        sta = np.load(allout)
        nameO,latO,lonO,durO = sta['NAME'],sta['LAT'],sta['LON'],sta['DURATION']

        # only certain event station pairs exist for ITO data
        if choice == "ITO":
            mseedpath = pathnames()['mseeds'] + 'ITO/{}*'.format(event_id)
            eventfiles = glob.glob(mseedpath)
            names = []
            for fid in eventfiles:
                names.append(fid.split('_')[1])

        # for ITO and LOBS data,. the process is the same from here except names
        newnames,newlats,newlons,newdurs = [],[],[],[]
        for i,sta in enumerate(names):
            try:
                if true_if_outside_bounds(lats[i],lons[i],corners):
                    print('\t {} is outside bounds'.format(sta))
                    continue
                
                ind = np.where(np.array(names)==sta)[0][0]
                if choice == "YH":
                    if 'EBS' in sta:
                        continue
                    code = "YH.{}.*.?H?".format(sta)
                else:
                    code = sta
                newnames.append(names[ind])
                newlats.append(lats[ind])
                newlons.append(lons[ind])
                duration = process_and_plot_waveforms(event_id,code,
                                                        choice=choice,
                                                        show=show,
                                                        save=save)
                newdurs.append(duration)

            except Exception as e:
                print('\t'+sta)
                newdurs.append(np.nan)
                plt.close()
                continue

        # output arrays are cats of old and new
        nameO = np.concatenate((nameO,np.array(newnames)))
        latO = np.concatenate((latO,np.array(newlats)))
        lonO = np.concatenate((lonO,np.array(newlons)))
        durO = np.concatenate((durO,np.array(newdurs)))

        dictout = {"NAME":nameO,"DURATION":durO,"LAT":latO,"LON":lonO}

    np.savez(allout,**dictout)


if __name__ == "__main__":
    event_id = '2014p715167'#'2014p864702
    global bounds
    bounds = [10,100]
    corners = [-42,-36,173,179.5]

    # for choice in ["YH"]:#["GN","ITO","YH"]:
    #     print(choice)
    #     loop_waveform_plotter(event_id,corners,choice,show=True,save=False)
    generate_duration_map(corners,event_id)
