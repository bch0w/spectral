"""test individual functions without having to run through the full script
"""
import sys
import numpy as np
sys.path.append('..')
import adjointBuilder
from obspy import read, read_events, read_inventory


def load_n_convert(fid):
    """windows saved as .npz are turned into numpy objects, rearrange structure
    so that they are the same as the output of run_pyflex
    """
    npzfile = np.load(fid)
    keys = npzfile.keys()
    output = {}
    for key in keys:
        output[key] = npzfile[key].tolist()
        
    return output
    
    
def test_calculate_adj_src():
    """testing creation of adjoint source using a multitaper approach
    """
    EXAMPLE_PAR_DICT = {"network":'NZ',
                        "station":'BFZ',
                        "code":'NZ.BFZ.*.HH?',
                        "event_id":'2014p240655',
                        "bounds":(6,30),
                        "rotate":True,
                        "units":'DISP',
                        "pyflex_config":'UAF',
                        "adj_src_type":'multitaper_misfit',
                        "save_adj_src":(False, None),
                        "comp_list":['R','T','Z'],
                        "save_plot":(False, None),
                        "plot":(True,False,True),
                        "dataset":None,
                        "verbose":True
                        }    
    st = read('./tests/test_data/testmseed.pickle')
    npzwindows = np.load('./tests/test_data/testwindows.npz')
    windows = convert_from_npz_window(npzwindows)
    
    adj_src = adjointBuilder.run_pyadjoint(st,windows,EXAMPLE_PAR_DICT)
    
                                 
def test_build_figure():
    """test figure building with example data
    """
    EXAMPLE_PAR_DICT = {"network":'NZ',
                        "station":'BFZ',
                        "code":'NZ.BFZ.*.HH?',
                        "event_id":'2014p240655',
                        "bounds":(6,30),
                        "rotate":True,
                        "units":'DISP',
                        "pyflex_config":'UAF',
                        "adj_src_type":'multitaper_misfit',
                        "save_adj_src":(False, None),
                        "comp_list":['R','T','Z'],
                        "save_plot":(False, None),
                        "plot":(True,False,True),
                        "dataset":None,
                        "verbose":True,
                        "stalta_wl":0.18
                        }    
    st = read('./tests/test_data/testmseed.pickle')
    cat = read_events('./tests/test_data/testevent.xml')
    event = cat[0]
    inv = read_inventory('./tests/test_data/testinv.xml')
    
    windows = load_n_convert('./tests/test_data/testwindows.npz')
    staltas = load_n_convert('./tests/test_data/teststalta.npz')
    adj_src = load_n_convert('./tests/test_data/testadjsrc.npz')

    adjointBuilder.build_figure(st,inv,event,windows,staltas,adj_src,
                                                            EXAMPLE_PAR_DICT)


def test_window_saving():
    """saving windows in pyasdf format as auxiliary data
    """
    windowpath = './tests/test_data/testwindows.npz'
    windows = np.load(windowpath)
    datasetname = pathnames()['data'] +'KUPEDATA/PYASDF/2014p240655.h5'
    DATASET = pyasdf.ASDFDataSet(datasetname,compression="gzip-3")
    for comp in windows.keys():
        internalpath = "{evid}/{net}_{sta}_{comp}".format(
                                                    evid='2014p240655',
                                                    net='NZ',
                                                    sta='BFZ',
                                                    comp=comp
                                                    )
        DATASET.add_auxiliary_data(data=windows[comp],
                                     data_type="PyflexWindows",
                                     path=internalpath,
                                     parameters={'test':'test'}
                                     )
                                     
