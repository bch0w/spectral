"""given a geonet event id, give all available station
"""

def accumulate_station_dates():
    """one-time script to grab start and endtimes for all stations to be used
    """
    from obspy import UTCDateTime, read_inventory
    from obspy.clients.fdsn import Client
    
    theDict = {}
    NAMES,STARTS,ENDS = [],[],[]
    
    # GEONET BROADBANDS
    c = Client("GEONET")
    new_zealand = [-50,-35,165,180]
    lat_lon = new_zealand
    inv = c.get_stations(network='NZ',
                        station='*Z',
                        channel='HH*',
                        minlatitude=lat_lon[0],
                        maxlatitude=lat_lon[1],
                        minlongitude=lat_lon[2],
                        maxlongitude=lat_lon[3],
                        level="station")

    # HOBITSS
    c = Client("IRIS")
    inv += c.get_stations(network='YH',
                        station="LOBS*",
                        location='',
                        level="station")
    inv += c.get_stations(network='YH',
                        station="EBS*",
                        location='',
                        level="station")  
                        
    # RDF
    inv += read_inventory(
                    '/Users/chowbr/Documents/subduction/RDF/XX.RDF.DATALESS')
    
    for network in inv:
        for station in inv[0]:
            NAMES.append("{n}.{s}".format(n=network.code,s=station.code))
            STARTS.append(station.start_date)
            ENDS.append(station.end_date)
    
    # SAHKE - broadbands not fdsn available, manual set times 
    SAHKE = {"LE4":[UTCDateTime('2010-136'),UTCDateTime('2010-331')],
            "LTN6":[UTCDateTime('2010-193'),UTCDateTime('2010-349')],
            "T004":[UTCDateTime('2010-088'),UTCDateTime('2010-255')],
            "T007":[UTCDateTime('2010-041'),UTCDateTime('2010-123')],
            "T010":[UTCDateTime('2010-135'),UTCDateTime('2010-348')],
            "T014":[UTCDateTime('2010-034'),UTCDateTime('2010-350')],
            "T016":[UTCDateTime('2010-088'),UTCDateTime('2010-322')],
            "T018":[UTCDateTime('2010-055'),UTCDateTime('2010-349')],
            "T020":[UTCDateTime('2010-089'),UTCDateTime('2010-261')] # poor avl
            }
            
    for station,timelist in SAHKE.items():
        NAMES.append(station)
        STARTS.append(timelist[0])
        ENDS.append(timelist[1])
    
    # BANNISTER - Information available from email correspondence
    BAN_path = ('/Users/chowbr/Documents/subduction/spectral/common/DATA/'
                'STATIONXML/Bannister_temporary_array_station_list.txt')
    with open(BAN_path,'r') as f:
        lines = f.readlines()
    
        
                        

    
    
    
    
      
    

        
        
    