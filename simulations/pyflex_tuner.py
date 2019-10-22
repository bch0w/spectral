"""
Tuning parameters in Pyflex is extremely important as bad measurements will lead
to bad kernels and therefore bad model updates.
"""
import os
import sys
import pyatoa
import pyasdf
import pyflex
import numpy as np
import logging

from scipy import signal


def set_pyflex_config(min_period, max_period, stalta_waterlevel=0.07, 
                      tshift_acceptance_level=10.,  dlna_acceptance_level=1.3, 
                      cc_acceptance_level=0.7,  s2n_limit=1.5, 
                      min_surface_wave_velocity=3., 
                      max_time_before_first_arrival=50., c_0=1., c_1=1.5, 
                      c_2=0., c_3a=4., c_3b=2.5, c_4a=2., c_4b=5., 
                      check_global_data_quality=False, 
                      snr_integrate_base=False, snr_max_base=3.):
    """
    """
    pyflex_config = pyflex.Config(
        min_period=min_period, max_period=max_period,
        stalta_waterlevel=stalta_waterlevel,
        tshift_acceptance_level=tshift_acceptance_level, tshift_reference=0.,
        dlna_acceptance_level=dlna_acceptance_level, dlna_reference=0.,
        cc_acceptance_level=cc_acceptance_level, s2n_limit=s2n_limit,
        earth_model="ak135",
        min_surface_wave_velocity=min_surface_wave_velocity,
        max_time_before_first_arrival=max_time_before_first_arrival,
        c_0=c_0, c_1=c_1, c_2=c_2, c_3a=c_3a, c_3b=c_3b, c_4a=c_4a, c_4b=c_4b,
        check_global_data_quality=check_global_data_quality,
        snr_integrate_base=snr_integrate_base, snr_max_base=snr_max_base,
        noise_start_index=0, signal_start_index=None, signal_end_index=-1,
        window_weight_fct=None,  window_signal_to_noise_type="amplitude",
        resolution_strategy="interval_scheduling"
        )

    return pyflex_config


def run_pyatoa(config, runtype="single", station="NZ.BFZ", save=False):
    """
    """
    # append title to title of figure
    append_title = ("\nstalta_waterlevel={0}; tshift_acceptance_level={1}; "
                    "dlna_acceptance_level={2}; cc_acceptance_level={3}; "
                    "\ns2n_limit={4}; max_time_before_first_arrival={5}; "
                    "c_0={6}; c_1={7}; c_2={8}; c_3a={9}; c_3b={10}; "
                    "\nc_4a={11}; c_4b={12}")
    pfc = config.pyflex_config[1]
    append_title = append_title.format(
            pfc.stalta_waterlevel, pfc.tshift_acceptance_level,
            pfc.dlna_acceptance_level, pfc.cc_acceptance_level,
            pfc.s2n_limit, pfc.max_time_before_first_arrival,
            pfc.c_0, pfc.c_1, pfc.c_2, pfc.c_3a, pfc.c_3b, 
            pfc.c_4a, pfc.c_4b
            )

    
    mgmt = pyatoa.Manager(config=config, empty=True)

    ds_path = os.path.join(os.getcwd(), "{}.h5".format(config.event_id))
    with pyasdf.ASDFDataSet(ds_path) as ds:
        mgmt.event = ds.events[0]
        if runtype == "single":
            sta_list = [station]
        else:
            sta_list = ds.waveforms.list()
        for sta in sta_list:
            mgmt.st_obs = ds.waveforms[sta].observed
            mgmt.st_syn = ds.waveforms[sta].synthetic_m00
            mgmt.inv = ds.waveforms[sta].StationXML
            mgmt.preprocess()
            mgmt.run_pyflex()
            mgmt.run_pyadjoint()
            mgmt.plot_wav(show=True, save=save, append_title=append_title, dpi=75)
            mgmt.reset(hard_reset=False)


def initialize(event_id, model_number):
    """
    """
    logger = logging.getLogger("pyflex")
    logger.setLevel(logging.DEBUG)

    config = pyatoa.Config(
                event_id=event_id,
                model_number=model_number,
                min_period=10,
                max_period=30,
                filter_corners=4,
                rotate_to_rtz=False,
                zero_pad=30,
                unit_output="DISP",
                window_amplitude_ratio=0.,
                pyflex_config="hikurang",
                adj_src_type="cc_hikurangi_strict",
                synthetics_only=True,
                cfgpaths={'synthetics':[],
                          'waveforms':[],
                          'responses':[]
                          }
                )

    figure_directory = os.path.join(os.getcwd(), model_number) 
    if not os.path.exists(figure_directory):
        os.makedirs(figure_directory)

    return config


def set_all_config():
    """
    """
    pyatoa_config = initialize(event_id="2013p142607", model_number="m00")
    pyflex_cfg = set_pyflex_config(
           min_period=pyatoa_config.min_period,
           max_period=pyatoa_config.max_period,
           stalta_waterlevel=0.07,
           s2n_limit=10.,
           tshift_acceptance_level=8.,
           dlna_acceptance_level=1.3,
           cc_acceptance_level=0.7,
           c_0=1.,
           c_1=1.5,
           c_2=0.,
           c_3a=4.,
           c_3b=2.5,
           c_4a=2.,
           c_4b=5.,
           min_surface_wave_velocity=3.,
           max_time_before_first_arrival=50.,
           check_global_data_quality=False,
           snr_integrate_base=3.5,
           snr_max_base=3.,
           )
    while True:
        while True:
            par = input("parameter= ")
            if not par:
                break
            val = input("value= ")
            setattr(pyflex_cfg, par, val)

        pyatoa_config.pyflex_config = ("test_config", pyflex_cfg)
        run_pyatoa(pyatoa_config, runtype="single", station="NZ.BFZ", save=False)

def iterate_parameters():
    """
    """
    pyatoa_config = initialize(event_id="2012p923684", model_number="m00")
    for param in np.arange(0.04, 0.1, 0.01):
        pyflex_cfg = set_pyflex_config(min_period=pyatoa_config.min_period,
                                       max_period=pyatoa_config.max_period,
                                       stalta_waterlevel=param
                                       )

        pyatoa_config.pyflex_config = ("test_config", pyflex_cfg)
        run_pyatoa(pyatoa_config, save=False)

        

if __name__ == "__main__":
    set_all_config()
    # iterate_parameters()


