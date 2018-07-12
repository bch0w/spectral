

def loop_itodata(event_id):
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
    pathout = pathnames()['data'] + 'DURATIONS'
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

def loop_lobsdata(event_id):
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
    pathout = pathnames()['data'] + 'DURATIONS'
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
