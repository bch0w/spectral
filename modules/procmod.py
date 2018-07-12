"""module file for processing data and data manipulation
"""
def myround(x,base=5,choice='near'):
    """round value x to nearest base, round 'up','down' or to 'near'est base
    """
    import numpy as np
    if choice == 'near':
        roundout = int(base * round(float(x)/base))
    elif choice == 'down':
        roundout = int(base * np.floor(float(x)/base))
    elif choice == 'up':
        roundout = int(base * np.ceil(float(x)/base))

    return roundout

def preprocess(st,resample=50,inv=None,output="VEL",filter=False):
    """preprocess waveform data:
    resample, demean, detrend, taper, remv. resp. (if applicable)
    """

    st_manipulate = st.copy()
    st_manipulate.resample(resample)
    st_manipulate.detrend("linear")
    st_manipulate.detrend("demean")
    st_manipulate.taper(max_percentage=0.05)
    inv_print = False
    if inv:
        inv_print = True
        st_manipulate.attach_response(inv)
        st_manipulate.remove_response(output=output,
                                      # pre_filt=[1/100,1/90,25,30],
                                      water_level=60,
                                      plot=False)
        st_manipulate.detrend("linear")
        st_manipulate.detrend("demean")
        st_manipulate.taper(max_percentage=0.05)

    # if no inventory, assumed to be synthetic data
    # change units from default velocity if necessary
    elif not inv:
        if output == "DISP":
            st_manipulate.differentiate(method="gradient")
        elif output == "ACC":
            st_manipulate.integrate(method="cumtrapz")
        st_manipulate.taper(max_percentage=0.05)

    code = st[0].get_id()
    print("\t[procmod.preprocess] {ID} {r}Hz resample, response: {i}".format(
                                                                ID=code,
                                                                r=resample,
                                                                i=inv_print))


    return st_manipulate

def signal_to_noise(data,separation):
    """calculate signal to noise ratio for a datastream
    """
    noise = data[:separation]
    signal = data[separation:]
    amplitude_SNR = abs(signal.max())/abs(noise.max())
    amplitude_SNR = round(amplitude_SNR,2)

    return amplitude_SNR

def trimstreams(st):
    """trim streams to common start and end times
    """
    st_trimmed = st.copy()
    start_set,end_set = 0,1E10
    for tr in st_trimmed:
        start_hold = tr.stats.starttime
        end_hold = tr.stats.endtime
        if start_hold > start_set:
            start_set = start_hold
        if end_hold < end_set:
            end_set = end_hold

    st_trimmed.trim(start_set,end_set)
    st_trimmed.detrend("linear")
    st_trimmed.detrend("demean")
    st_trimmed.taper(max_percentage=0.05)

    return st_trimmed

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



if __name__ == "__main__":
    print('what?')
