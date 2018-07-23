"""misfit visualization tool to be called through adjointBuilder
produces a basemap with beachball and all available stations as well as the
relevant station highlighted. important information annotated (e.g. 
misift information, distance, BAz etc.)
"""

def event_beachball(MT,shift=[0,0]):
    """plot event beachball on basemap 'm' object for a given geonet event_id
    """
    from obspy.imaging.beachball import beach
    eventx,eventy = m(MT['Longitude'],MT['Latitude'])
    FM = [MT['strike2'],MT['dip2'],MT['rake2']]

    b = beach(FM,xy=(eventx,eventy),
                 width=2.5E4,
                 linewidth=1,
                 facecolor='r')
    b.set_zorder(1000)
    ax = plt.gca()
    ax.add_collection(b)
    yshift = shift[0] * (m.ymax-m.ymin)
    xshift = shift[1] * (m.xmax-m.xmin)
    plt.annotate("M{m}\n{e}\n(depth: {d:} km)".format(m=MT['Mw'],
                                         e=event_id,#!!!
                                         d=int(MT['CD'])),
                 xy=(eventx,eventy),
                 xytext=(eventx+xshift,eventy+yshift),
                 fontsize=10,
                 zorder=200,
                 weight='bold',
                 multialignment='center')


def generate_duration_map(corners):
    """initiate and populate a basemap object for New Zealands north island.
    Functionality to manually ignore stations based on user quality control
    Takes station coordinates and coloring from npz files
    Choice to annotate two stations which correspond to waveforms
    Calls beachball and trench tracer to populate basemap
    :type corners: list of floats
    :param corners: values for map corners to set bounds
     e.g. [lat_bot,lat_top,lon_left,lon_right]
    """
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
        if MI in names:
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
            continue
            color = 'w'
            D = 10
        else:
            color = colormap.to_rgba(durs[i])
        X,Y = m(lons[i],lats[i])
        m.scatter(X,Y,marker='v',
                        color=color,
                        s=100,
                        linewidth=1.1,
                        edgecolor='k',
                        zorder=int(D))

        # fine tuned annotations
        """annotation dictionary
        2014p051675 HAZ yshift=.0035 xshift=.0125, PUZ yshift=-.0255 xshift=-.05
        2015p768477 TOZ yshift=.01 xshift=.01, KU15-3 yshift=-.0255 xshift=-.075
        2014p715167 OPRZ yshift=.01 xshift=.01 , PUZ yshift=.015 xshift=-.05
        2014p240655 MXZ yshift=.01 xshift=.01, MKAZ yshift=-.025,xshift=-.05
        2014p864702 MXZ yshift=.01 xshift=.015, EBPR-2 yshift=.005, xshift=.015
        """
        ano = {"2014p051675":["HAZ",.0035,.0125,"PUZ",-.0255,-.05],
               "2015p768477":["TOZ",.01,.01,"KU15-3",-.0255,-.075],
               "2014p864702":["MXZ",.01,.015,"EBPR-2",.005,.015],
               "2015p822263":["MRZ",.01,.015,"PUZ",-.0255,0.015],
               "2014p240655":["MXZ",.01,.015,"MKAZ",-.025,-.05]
               }
               
        STA1 = ano[event_id][0]
        STA2 = ano[event_id][3]
        if (names[i] == STA1) or (names[i] == STA2):
            if names[i] == STA1:
                yshift = ano[event_id][1] * (m.ymax-m.ymin)
                xshift = ano[event_id][2] * (m.xmax-m.xmin)
            elif names[i] == STA2:
                yshift = ano[event_id][4] * (m.ymax-m.ymin)
                xshift = ano[event_id][5] * (m.xmax-m.xmin)
            plt.annotate(names[i],
                         xy=(X,Y),
                         xytext=(X+xshift,Y+yshift),
                         fontsize=10,
                         zorder=1000,
                         weight='bold')
        # plt.annotate(names[i],
        #              xy=(X,Y),
        #              xytext=(X,Y),
        #              fontsize=10,
        #              zorder=1000,
        #              weight='bold')

    plt.annotate("(depth: 21 km)\n(depth: 14 km)\n(depth: 20 km)\n(depth: 60 km)\n(depth: 19 km)",
                xy=(X,Y+0.02*(m.ymax-m.ymin)),fontsize=8,zorder=1000,weight='bold')
    # additional map objects
    trace_trench(m)
    event_beachball_durationmap(event_id,m,anno=True)
    # m.drawmapscale(179.65,-41.75, 179.65,-41.75, 100,yoffset=0.01*(m.ymax-m.ymin))
    m.drawmapscale(179,-41.75, 179,-41.75, 100,yoffset=0.01*(m.ymax-m.ymin))

    # colorbar
    colormap.set_array(durs)
    cbar = f.colorbar(colormap,fraction=0.046,pad=0.04)
    cbar.set_label('duration (s)',rotation=270,labelpad=17,fontsize=14.5)

    # plt.title(event_id)

    if save:
        outputfolder = pathnames()['spectralplots'] + 'durations/PAPEROUT'
        fid = 'map_{e}_{b0}_{b1}.png'.format(e=event_id,
                                                b0=bounds[0],b1=bounds[1])
        fidout = os.path.join(outputfolder,fid)
        plt.savefig(fidout,dpi=plt.gcf().dpi)
    if show:
        plt.show()