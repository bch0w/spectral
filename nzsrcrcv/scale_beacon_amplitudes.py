def scale_beacon_amplitudes(st, st_syn):
    """
    Warning:
        Not actually used, we have decided that the affected intruments should
        just be ignored rather than hacked.

    Return a waveform with scaled amplitudes based on station identifiers.
    This scaling is empirical and is required to boost amplitudes of our
    CMG-40T60s instruments which for some reason have order of magnitude lower
    amplitudes. The 30s instruments are okay but amplitudes are larger than
    GeoNet stations (maybe site response?), so we put in a small scaling factor.

    Station Number for 60s instruments (require scaling):
        1, 2, 4, 5,  6, 7, 8, 9, 13, 14, 15, 18
    Station Number for 30s instruments:
        3, 10, 11, 12, 16, 17, 19, 20, 21
    :type st: obspy.core.stream.Stream
    :param st: stream containing waveform data to be scaled
    :type st_syn: obspy.core.stream.Stream
    :param st_syn: synthetic waveform used to check the scaled amplitudes
    :rtype: obspy.core.stream.Stream
    :return: stream with scaled waveforms
    """
    st_scale = st.copy()

    # Amplitude scaling
    scale_60s = 75.
    instr_60s = [1, 2, 4, 5, 6, 7, 8, 9, 13, 14, 15, 18]

    # Determine group by checking station name
    try:
        idx = int(st_scale[0].stats.station[2:])
    except ValueError:
        logger.debug("station code does not match BEACON station formatting")
        return st

    # Scale the group based on the instrument type
    if idx in instr_60s:
        logger.debug(f"scaling {st_scale[0].get_id()} by -1 * {scale_60s}")
        for tr in st_scale:
            data_try = tr.data * scale_60s
            # Get the synthetic component to use for comparisons
            comp = tr.id[-1]
            tr_syn = st_syn.select(component=comp)[0]
            # Check the scaled data
            if data_try.max() > 10 * tr_syn.data.max():
                logger.warning(f"scaled obs max is >10x syn max")
            # Check if the scaled data
            tr.data = data_try

    return st_scale

