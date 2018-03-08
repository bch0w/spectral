import sys
sys.path.append('../modules/')
from getdata import pathnames
from synmod import stf_convolve, preprocess

def preprocess(st,inv):
    """preprocess waveform data
    """
    st_manipulate = st.copy()
    st_manipulate.detrend("demean")
    st_manipulate.detrend("linear")
    st_manipulate.taper(max_percentage=0.05)
    st_manipulate.attach_response(inv)
    st_manipulate.remove_response(output="VEL",water_level=0)

    return st_manipulate
