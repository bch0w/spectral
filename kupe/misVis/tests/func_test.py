"""test individual functions without having to run through the full script
"""
import sys
import numpy as np
sys.path.append('..')
import adjointBuilder
from obspy import read, read_events, read_inventory


def convert_from_npz_window(npzwindow):
    """windows saved as .npz are turned into numpy objects, rearrange structure
    so that they are the same as the output of run_pyflex
    """
    keys = npzwindow.keys()
    windows = {}
    for key in keys:
        windows[key] = npzwindow[key].tolist()
        
    return windows
    
    
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
    boundsdict = {"station_name":"TEST","bounds":(6,30)}
    streampath = './tests/test_data/testmseed.pickle'
    windowpath = './tests/test_data/testwindows.npz'
    st = read(streampath)
    windows = np.load(windowpath)
    eventpath = './tests/test_data/testevent.xml'
    invpath = './tests/test_data/testinv.xml'
    cat = read_events(eventpath)
    event = cat[0]
    inv = read_inventory(invpath)


    build_figure(st,inv,event,windows,boundsdict)


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
                                     
