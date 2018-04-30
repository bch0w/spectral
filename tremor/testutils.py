"""functions used to test outputs and code segments of pyfreqscan
"""
def test_TEROR_creation(st_raw, inv, night, already_preprocessed=False):
    """create passband arrays more in line with methods by Sit et al. 2012
    :type st_raw: obspy stream
    :param st_raw: raw stream
    :type inv: obspy inventory
    :param inv: response information for st_raw
    :rtpye TEORRm: list of numpy arrays
    :return TEORRm: arrays containing filtered waveforms and amp. ratios
    """
    # parameter set
    tremor_band = [2,8]
    tremor_band_short = [2,5]
    earthquake_band = [10,20]
    earthquake_band_short = [10,15]
    ocean_band = [0.5,1]
    surface_wave_band = [0.02,0.1]

    sampling_rate = 50
    five_seconds = sampling_rate * 5
    one_minute = sampling_rate * 60
    five_minutes = one_minute * 5
    one_hour = one_minute * 60
    
    # preprocess
    if not already_preprocessed:
        st = preprocess(st_raw,inv,resample=sampling_rate,night=night)
    else:
        st = st_raw
    # create arrays for different freq bands
    tremor_horizontal = create_horizontal_data(st,tremor_band)
    tremor_horizontal_short = create_horizontal_data(st,tremor_band_short)
    earthquake_horizontal = create_horizontal_data(st,earthquake_band)
    earthquake_horizontal_short = create_horizontal_data(st,
                                                    earthquake_band_short)
    surface_wave_horizontal = create_horizontal_data(st,surface_wave_band)
    ocean_horizontal, water_level = set_water_level(st,ocean_band)
    
    # create normalized surface_wave amplitudes
    S = average_array(surface_wave_horizontal,one_hour)
    S_min = surface_wave_horizontal.min()
    S /= S_min
    

    
    import matplotlib.pyplot as plt
    f,(ax1,ax2,ax3,ax4) = plt.subplots(4,sharex=False,sharey=False,
                                        figsize=(11.69,8.27),dpi=75)
    ax1.plot(tremor_horizontal,'k')
    ax1.plot(tremor_horizontal_short,'r')
    ax2.plot(earthquake_horizontal,'k')
    ax2.plot(earthquake_horizontal_short,'r')
    ax3.plot(surface_wave_horizontal,'k')
    ax3.plot(ocean_horizontal,'r',zorder=2)
    ax4.plot(S,'k')
    plt.show()
    sys.exit()
    
    
    # create ratio arrays by certain time windows
    Rs,Rm,Rh = [],[],[]
    ratio_equation = lambda T,E,O: T**2 / (E*O)
    for window,Rx in zip([five_seconds,five_minutes,one_hour],[Rs,Rm,Rh]):
        T = average_array(tremor_horizontal,window)
        E = average_array(earthquake_horizontal,window)
        S = average_array(surface_wave_horizontal,window)
        R = ratio_equation(T=T,E=E,O=S)
        Rx.append(R)