"""module file for processing data
"""
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
<<<<<<< HEAD
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
        
=======
        st_manipulate.taper(max_percentage=0.05)

>>>>>>> 88adbac2105e0c3c1e91e5c44549b40695098512
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

def amplitude_threshold(t,tr,threshold_percentage):
    """based on some amplitude threshold, return all points that fall above,
    and the time length covered
    """
    duration_a,tover_plot,aover_plot,threshold_plot,sample_plot = [],[],[],[],[]

    samp_rate = tr.stats.sampling_rate
    tr_manipulate = (tr.data**2)**(1/2)
    peak_amp = tr_manipulate.max()
    threshold = peak_amp * threshold_percentage
    threshold_plot.append(threshold)

    # loop over seismogram, determine start and end of peak energy
    # a for amplitude, s for sample
    a_over, s_over = [],[]
    for i,amplitude in enumerate(tr_manipulate):
        if amplitude >= threshold:
            s_over.append(i)
            a_over.append(amplitude)

    # find edgepoints by checking if the next sample j is the same as i+1
    s_edge,a_edge,sections,samples = [s_over[0]],[a_over[0]],[],[]
    for i,(S,A) in enumerate(zip(s_over[1:-2],a_over[1:-2])):
        if s_over[i+2] != (S + 1):
            # determine number of samples covered
            samples.append(len(s_edge))
            s_edge,a_edge = [s_over[i+1]],[a_over[i+1]]
        else:
            # if the next sample is the same, keep going
            s_edge.append(S)
            a_edge.append(A)

    duration_in_seconds = round(sum(samples)/samp_rate,0)

    # convert samples to time
    t_over = []
    for S in s_over:
        t_over.append(t[S])

    return t_over,a_over,duration_in_seconds


if __name__ == "__main__":
    print('what?')
