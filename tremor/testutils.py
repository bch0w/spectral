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
        

def tremor_counter(tc_array,choice,avg_choice='mean',night=False):
    """count number of tremors per day using standard deviation thresholds
    :type tc_array: np.array
    :param tc_array: average values for amplitude ratios
    :type nighttime: bool
    :param nighttime: only nighttime, files with '_night' suffix
    :rtype sigma_dict: dict of np.arrays
    :return sigma_dict: 2 and 3-sigma arrays
    """
    sig_dict = mean_std_creator(night)
    mean_ = sig_dict["{}_mean".format(choice)]
    median_ = sig_dict["{}_median".format(choice)]
    one_sigma = sig_dict["{}_sigma".format(choice)]

    avg_val = mean_
    if avg_choice == 'median':
        avg_val = median_

    two_sigma = avg_val + one_sigma * 2
    three_sigma = avg_val + one_sigma * 3

    R_2sig = np.copy(tc_array)
    R_3sig = np.copy(tc_array)

    # count tremor activity by fulfilling threshold criteria for 2sig and 3sig
    for i,section in enumerate(tc_array):
        if section < two_sigma:
            R_2sig[i] = 0
            R_3sig[i] = 0
        elif section < three_sigma:
            R_3sig[i] = 0

    print("[tremor_counter {c}: 2-sig:{two} / 3-sig:{three}]".format(c=choice,
                                                    two=len(R_2sig[R_2sig>0]),
                                                    three=len(R_3sig[R_3sig>0])
                                                    ))

    sigma_dict = {"2-sigma":R_2sig,
                  "3-sigma":R_3sig}

    return sigma_dict

def count_tremors_over_sigma(TEORRm,avg_choice="median",night=False,
                                               stop_if_tremor_num_below = 3):
    """count number of tremors in a given array and return values which cross
    a global 2-sigma and 3-sigma threshold, which signifies a tremor 'detection'
    station specific function. Currently only returns 2-sigma values
    :return tremor_count_dict: dictionary filled with arrays where ratio arrays
    surpass 2-sigma over the mean value
    """
    tremor_count_dict = {}
    for choice in ['Rs','Rm','Rh']:
        sgtmp = tremor_counter(TEORRm[choice],choice=choice,
                                        avg_choice=avg_choice,night=night)
        two_sigma = sgtmp["2-sigma"]
        ^R_2sig
        # check-stop if low-tremor detection threshold met on Rs
        if stop_if_tremor_num_below and choice=='Rm':
            if len(two_sigma[two_sigma>0]) < stop_if_tremor_num_below:
                print("\t--tremor number below {n} ".format(
                                                n=stop_if_tremor_num_below))
                return None

        # nan out 0's to avoid overplotting 0 values
        two_sigma[two_sigma==0] = np.nan
        tremor_count_dict[choice] = two_sigma

    return tremor_count_dict
