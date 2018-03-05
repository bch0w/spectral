"""module containing functions pertaining to synthetic outputs of specfem
"""
def mt_transform(mt,method):
    """transform moment tensor between xyz and rtp coordinates
    acceptable mt formats:
        [m11,m22,m33,m12,m13,m23]
        [mxx,myy,mzz,mxy,mxz,myz]
        [mrr,mtt,mpp,mrt,mrp,mtp]
    :type mt: list
    :param mt: moment tensor in format above
    :type method: str
    :param method: type of conversion, "rtp2xyz" or "xyz2rtp"
    """
    if method == 'xyz2rtp':
        m_rr = mt[2] #m_zz
        m_tt = mt[0] #m_xx
        m_pp = mt[1] #m_yy
        m_rt = mt[4] #m_xz
        m_rp = -mt[5] #m_yz
        m_tp = -mt[3] #m_xy
        return [m_rr,m_tt,m_pp,m_rt,m_rp,m_tp]
    if method == 'rtp2xyz':
        m_xx = mt[1] #m_tt
        m_yy = mt[2] #m_pp
        m_zz = mt[0] #m_rr
        m_xy = -mt[5] #m_tp
        m_xz = mt[3] #m_rt
        m_yz = -mt[4] #m_rp
        return [m_xx,m_yy,m_zz,m_xy,m_xz,m_yz]
    else:
        print("Invalid transformation method")

def generate_CMTSOLUTION(event_id):
    """generate CMTSOLUTION file in the format of the Harvard CMT catalog
    using values from Ristau's moment tensor solutions
    """
    import sys
    sys.path.append("../spectral/")
    from obspy import UTCDateTime
    from getdata import get_moment_tensor
    from obspy.clients.fdsn import Client
    from obspy.geodetics import FlinnEngdahl

    c = Client('IRIS')

    # moment tensor dictionary, put tensor in proper units, prepare entries
    MT = get_moment_tensor(event_id=event_id)
    mt = [MT['Mxx'],MT['Myy'],MT['Mzz'],MT['Mxy'],MT['Mxz'],MT['Myz']]
    mt = [_*(1E20) for _ in mt]
    mt = mt_transform(mt,method='xyz2rtp')
    fe = FlinnEngdahl()
    region = fe.get_region(MT['Longitude'],MT['Latitude'])
    date = UTCDateTime(str(int(MT['Date'])))
    # first line
    header = "XXXX {Y} {M} {D} {H} {m} {S} {La} {Lo} {Dp} {Mb} {Ms} {N}".format(
            Y=date.year,
            M=date.month,
            D=date.day,
            H=date.hour,
            m=date.minute,
            S=date.second,
            La=MT['Latitude'],
            Lo=MT['Longitude'],
            Dp=MT['CD'],
            Mb=0,
            Ms=0,
            N=region)
    # write to solution file
    with open('./CMTSOLUTIONS/{}CMTSOLUTION'.format(event_id),'w') as f:
        f.write(header+'\n')
        f.write("event name:\t {}\n".format(event_id))
        f.write("time shift:\t 0\nhalf duration: 0\n")
        f.write("latorUTM:\t {}\n".format(MT['Latitude']))
        f.write("longorUTM:\t {}\n".format(MT['Longitude']))
        f.write("depth:\t{}\n".format(MT['CD']))
        f.write("Mrr:\t{}\n".format(mt[0]))
        f.write("Mtt:\t{}\n".format(mt[1]))
        f.write("Mpp:\t{}\n".format(mt[2]))
        f.write("Mrt:\t{}\n".format(mt[3]))
        f.write("Mrp:\t{}\n".format(mt[4]))
        f.write("Mtp:\t{}\n".format(mt[5]))



def stf_convolve(st,half_duration,window="bartlett"):
    """convolve source time function with a stream
    :type st: obspy.stream
    :param st: stream object containing traces of data
    :type half_duration: float
    :param half_duration: half duration of stf in seconds
    :type window: str
    :param window: window type to return
    ========================================
    boxcar, triang, blackman, hamming, hann,
    bartlett, flattop, parzen, bohman,
    blackmanharris, nuttall, barthann,
    kaiser (needs beta),
    gaussian (needs standard deviation),
    general_gaussian (needs power, width),
    slepian (needs width),
    chebwin (needs attenuation),
    exponential (needs decay scale),
    tukey (needs taper fraction)
    NOTE: bartlett window is a triangle that touches 0
    ========================================
    :return st:
    """
    import numpy as np
    from scipy import signal

    # set up window
    npts = st[0].stats.npts
    sampling_rate = st[0].stats.sampling_rate
    half_duration_in_samples = round(half_duration * sampling_rate)
    stf = signal.get_window(window=window,
                            Nx=half_duration_in_samples * 2)

    # convolve
    new_st = st.copy()
    for tr in new_st:
        new_data = np.convolve(tr.data,stf,mode="same")
        tr.data = new_data

    return new_st
