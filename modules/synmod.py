"""module containing functions pertaining to synthetic outputs of specfem
"""
import sys
import getdata
from getdata import pathnames
from obspy import read_events

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


def mt_from_event(event):
    """pull out the tensor components from obspy event object, return as dict
    """
    focmec = event.focal_mechanisms[0].moment_tensor.tensor
    m_rr = focmec.m_rr
    m_tt = focmec.m_tt
    m_pp = focmec.m_pp
    m_rt = focmec.m_rt
    m_rp = focmec.m_rp
    m_tp = focmec.m_tp

    MT = {"Mrr":m_rr,"Mtt":m_tt,"Mpp":m_pp,
          "Mrt":m_rt,"Mrp":m_rp,"Mtp":m_tp}

    return MT

def compare_beachballs(event_id):
    """plot beachballs to compare GCMT and converted GEONET MT
    """
    from obspy.imaging.beachball import beachball
    from getdata import get_moment_tensor
    import matplotlib.pyplot as plt

    MT = get_moment_tensor(event_id=event_id)
    geonet_xyz = [MT['Mxx'],MT['Myy'],MT['Mzz'],MT['Mxy'],MT['Mxz'],MT['Myz']]
    geonet_rtp = mt_transform(geonet_xyz,method='xyz2rtp')

    gcmt = getdata.get_GCMT_solution(event_id=event_id)
    GCMT = mt_from_event(gcmt)
    gcmt_rtp = [GCMT['Mrr'],GCMT['Mtt'],GCMT['Mpp'],
                GCMT['Mrt'],GCMT['Mrp'],GCMT['Mtp']]
    gcmt_xyz = mt_transform(geonet_rtp,method='rtp2xyz')

    bp = pathnames()['kupeplots'] + 'beachballs/{e}_{o}_{c}.png'
    
    geonet_xyz_b = beachball(geonet_xyz,outfile=bp.format(e=event_id,
                                                          o='geonet',
                                                          c='xyz',
                                                          ))
    geonet_rtp_b = beachball(geonet_rtp,outfile=bp.format(e=event_id,
                                                          o='geonet',
                                                          c='rtp',
                                                          ))
    gcmt_rtp_b = beachball(gcmt_rtp,outfile=bp.format(e=event_id,
                                                      o='gcmt',
                                                      c='rtp',
                                                      ))
    gcmt_xyz_b = beachball(gcmt_xyz,outfile=bp.format(e=event_id,
                                                      o='gcmt',
                                                      c='xyz',
                                                      ))
    plt.close("all")


def stf_convolve(st,half_duration,window="bartlett",time_shift=False):
    """convolve source time function with a stream, time shift if needed
    :type st: obspy.stream
    :param st: stream object containing traces of data
    :type half_duration: float
    :param half_duration: half duration of stf in seconds
    :type window: str
    :param window: window type to return
    :type time_shift: float
    :param time_shift: change the starttime
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
    :return new_st:
    """
    import numpy as np
    from scipy import signal

    # set up windowtriang
    npts = st[0].stats.npts
    sampling_rate = st[0].stats.sampling_rate
    half_duration_in_samples = round(half_duration * sampling_rate)
    stf = signal.get_window(window=window,
                            Nx=(half_duration_in_samples * 2) -1)

    # make sure window touches 0 at the end
    if stf[-1] != 0:
        stf = np.append(stf,0)

    # normalize to keep area of window equal one
    stf *= (2/len(stf))

    # convolve and time shift if necessary
    new_st = st.copy()
    for tr in new_st:
        if time_shift:
            tr.stats.starttime = tr.stats.starttime + time_shift
        new_data = np.convolve(tr.data,stf,mode="same")
        tr.data = new_data


    return new_st

def tshift_halfdur(event_id):
    """get the absolute time shift between centroid time and hypocenter time to
    shift the synthetic seismogram into absolute time using the equation
    t_abs = t_pde + time shift + t_syn
    also get the half duration from the GCMT solution
    """
    # time shift
    MT = getdata.get_GCMT_solution(event_id)
    CMTSOLUTIONPATH = (pathnames()['kupedata'] +
                                'CMTSOLUTIONS/{}CMTSOLUTION'.format(event_id))
    CMTSOLUTION = read_events(CMTSOLUTIONPATH)

    CMTSOLUTION_time = CMTSOLUTION[0].origins[0].time
    CENTROID_time = [i.time for i in MT.origins
                                            if i.origin_type == "centroid"][0]

    time_shift = abs(CMTSOLUTION_time - CENTROID_time)

    # half duration
    moment_tensor = MT.focal_mechanisms[0].moment_tensor
    half_duration = (moment_tensor.source_time_function['duration'])/2

    return time_shift, half_duration

if __name__ == "__main__":
    print('what?')
