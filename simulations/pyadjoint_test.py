"""logging.basicConfig(filename=logname,
                            filemode='a',
                            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.DEBUG)
Testing Pyflex configurations to determine best fitting parameters for optimal
window picking
"""
import pyasdf
import pyatoa
import logging
import traceback


def log(switch):
    """
    set logging on or off
    :param switch:
    :return:
    """
    if switch:
        logger_pyatoa = logging.getLogger("pyatoa")
        logger_pyatoa.setLevel(logging.DEBUG)
        
        logger_pyflex = logging.getLogger("pyflex")
        logger_pyflex.setLevel(logging.DEBUG)

        logger_pyadjoint = logging.getLogger("pyadjoint")
        logger_pyadjoint.setLevel(logging.DEBUG)


def short_run():
    """
    Quickly run through Pyatoa framework
    :return:
    """
    pyflex_config = "default"
    pyadjoint_branch = "dev"
    folder = "eightevent_m03s00"

    config = pyatoa.Config(
        event_id="2016p275188",
        model_number="m02",
        min_period=10,
        max_period=30,
        filter_corners=4,
        rotate_to_rtz=False,
        unit_output="DISP",
        pyflex_config=pyflex_config,
        adj_src_type="mtm_hikurangi_strict",
        synthetics_only=False,
        window_amplitude_ratio=0.2
        )

    mgmt = pyatoa.Manager(config=config, empty=True)

    with pyasdf.ASDFDataSet(
            "./{}/{}.h5".format(folder, config.event_id)) as ds:
        mgmt.event = ds.events[0]
        for sta in ds.waveforms.list():
            try:
                if True:
                # if sta in ["NZ.BKZ"]:
                    print(sta)
                    mgmt.st_obs = ds.waveforms[sta].observed
                    mgmt.st_syn = ds.waveforms[sta].synthetic_m00
                    mgmt.inv = ds.waveforms[sta].StationXML
                    mgmt.preprocess()
                    mgmt.run_pyflex()
                    mgmt.run_pyadjoint()
                    # mgmt.plot_map(show="hold")
                    mgmt.plot_wav(show=True,
                                  save="./figures/{}_{}_{}.png".format(
                                      sta, pyflex_config, pyadjoint_branch),
                                  append_title=pyflex_config
                                  )
                    mgmt.reset(hard_reset=False)
            except Exception as e:
                traceback.print_exc()
                import ipdb;ipdb.set_trace()
                print("error")
                mgmt.reset(hard_reset=False)
                continue


if __name__ == "__main__":
    log(True)
    short_run()
