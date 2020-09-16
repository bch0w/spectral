"""
Functionality to test windowing parameters using Pyatoa and Pyflex
"""
from pyasdf import ASDFDataSet
from pyatoa import Manager
from pyatoa.core.config import set_pyflex_config


class Info(dict):
    """
    Dictionary object that allows attribute selection and contains properties
    that relate to the various Pyflex windowing parameters.
    """
    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]

    def __repr__(self):
        """
        List the available properties by difference with inhereted properties
        :return:
        """
        properties = list(set(dir(self)) - set(dir(dict)))
        # Get rid of internal attributes
        properties = sorted([_ for _ in properties if not _.startswith("__")])
        return str(properties)

    @property
    def tshift_acceptance_level(self):
        print("""TIME SHIFT ACCEPTANCE [s]
                 Maximum allowable cross correlation time shift/lag relative to 
                 the reference. 
        """)
        return

    @property
    def dlna_acceptance_level(self):
        print("""AMPLITUDE RATIO ACCEPTANCE
                 Maximum allowable amplitude ratio
        """)
        return

    @property
    def cc_acceptance_level(self):
        print("""AMPLITUDE RATIO ACCEPTANCE
                 Maximum allowable normalized cross correlation per window
        """)
        return

    @property
    def c_0(self):
        print("""FOR REJECTION OF INTERNAL MINIMA
                 A multiple of water level wE. 
                 No window may contain a local 
                 minima in its STA:LTA waveform that falls below 
                 the local value of: c0 * wE
                 + Larger c_0: More conservative
        """)
        return

    @property
    def c_1(self):
        print("""FOR REJECTION OF SHORT WINDOWS
                 A multiple of min period T0. 
                 No window shorter than C1 * T0
                 + Larger c_1: Windows must be longer
        """)
        return

    @property
    def c_2(self):
        print("""FOR REJECTION OF UN-PROMINENT WINDOWS
                 A multiple of of water level wE. 
                 A window whose seed maximum on the STA:LTA waveform rises less
                 than C2 * wE above either of its adjacent minima is rejected.
                 NOTE: C2 is hard to control. Best left as 0.
        """)
        return

    @property
    def c_3a(self):
        print("""FOR REJECTION OF MULTIPLE DISTINCE ARRIVALS
                 Expressed as a fraction. Regulates acceptable height ratio 
                 between local maxima in a given window on the STA:LTA waveform.
                 + 
        """)
        return

    @property
    def c_3b(self):
        print("""FOR REJECTION OF MULTIPLE DISTINCE ARRIVALS
                 Expressed as a multiple of min period T0. Regulates acceptable 
                 height ratio between local maxima in a given window on the 
                 STA:LTA waveform.
        """)
        return

    @property
    def c4a(self):
        print("""FOR REJECTION OF EMERGENT STARTS AND/OR CODAS
                 Expressed as a multiple of min period T0. Limits the length
                 of a window before its first local maximum in STA:LTA
        """)
        return

    @property
    def c_4b(self):
        print("""FOR REJECTION OF EMERGENT STARTS AND/OR CODAS
                 Expressed as a multiple of min period T0. Limits the length
                 of a window after its last local maximum in STA:LTA
        """)
        return

    @property
    def s2n_limit(self):
        print("""SIGNAL-TO-NOISE LIMIT
                 Limit of the signal to noise ratio per window. If the maximum 
                 amplitude of the window over the maximum amplitude of the 
                 global noise of the waveforms is smaller than this window, 
                 then it will be rejected.
        """)
        return

    @property
    def check_global_data_quality(self):
        print("""CHECK GLOBAL DATA QUALITY
                 Determines whether or not to check the signal to noise ratio of 
                 the whole observed waveform. If True, no windows will be 
                 selected if the signal to noise ratio is above the thresholds.
        """)
        return

    @property
    def snr_integrate_base(self):
        print("""SIGNAL-TO-NOISE MINIMUM RATIO
                 Minimal SNR ratio. If the squared sum of the signal normalized 
                 by its length over the squared sum of the noise normalized by 
                 its length is smaller then this value, no windows will be 
                 chosen for the waveforms. 
                 Only used if check_global_data_quality is True.
        """)
        return

    @property
    def snr_max_base(self):
        print("""SIGNAL-TO-NOISE MAXIMUM RATIO
                 Minimal amplitude SNR ratio. If the maximum amplitude of the 
                 signal over the maximum amplitude of the noise is smaller than 
                 this value no windows will be chosen for the waveforms. 
                 Only used if check_global_data_quality is True.
        """)
        return


class WasherError(Exception):
    """Generic error here"""
    pass


class WindowWasher:
    """
    Test time windowing using Pyflex presets
    """
    def __init__(self, ds_fid):
        """
        The class contains a PyASDF Dataset that will be used to load observed
        and synthetic data, based on user input
        """
        self.ds = ASDFDataSet(ds_fid)
        # Initiate an empty Manager to get access to its config
        self.mgmt = Manager()
        self.info = Info()

    def setup(self, code):
        """
        Load waveforms from the dataset, and preprocess for windowing
        """
        # If Configure was run first, this will pass its config object over
        try:
            self.mgmt = Manager(ds=self.ds, config=self.mgmt.config)
            self.mgmt.load(code, config=False, synthetic_tag="synthetic_m00s00")
            self.mgmt.config.save_to_ds = False
            self.mgmt.standardize().preprocess()
        except Exception as e:
            raise WasherError(e)

    def configure(self, min_period, max_period, preset="default", **kwargs):
        """
        Configure the Pyflex Parameters, either by passing kwargs or
        choosing preset
        :return:
        """
        self.mgmt.config.min_period = min_period
        self.mgmt.config.max_period = max_period
        self.mgmt.config.pyflex_preset = preset
        self.mgmt.config.pyflex_config, unused = set_pyflex_config(
            choice=preset, min_period=min_period, max_period=max_period,
            **kwargs
        )
        if unused:
            for uu in unused:
                print(f"{uu} is not a valid Pyflex parameter")

    def set(self, **kwargs):
        """
        After a Pyflex Config has been set by 'configure', can 'set' individual
        parameters for fine-tuning
        """
        assert self.mgmt is not None
        for key in kwargs.keys():
            assert hasattr(self.mgmt.config.pyflex_config, key), f"{key} no no"

        vars(self.mgmt.config.pyflex_config).update(**kwargs)

    def check(self, key):
        """
        Check what the current kwarg value of the Pyflex config is
        """
        assert self.mgmt is not None
        if hasattr(self.mgmt.config.pyflex_config, key):
            print(f"{key}: {getattr(self.mgmt.config.pyflex_config, key)}")
        else:
            print(f"{key}: NOT FOUND")

    def wash(self):
        """
        Run Manager window() function and plot the waveforms, no map, to look at
        window quality.
        """
        assert(self.mgmt is not None)
        self.mgmt.window()
        self.mgmt.plot(choice="wav", show=True, save=None)


if __name__ == "__main__":
    # Call window washer
    import time
    from glob import glob
    from pyflex import logger
    from pyatoa.core.manager import ManagerError

    logger.setLevel("DEBUG")

    total = 0
    for fid in glob("*.h5"):
        ww = WindowWasher(ds_fid=fid)
        ww.configure(15, 30, "nznorth_15-30s")
        for sta in ww.ds.waveforms.list():
            try:
                ww.setup(sta)
            except WasherError:
                continue
            try:
                ww.wash()
                time.sleep(3)
                total += ww.mgmt.stats.nwin
            except ManagerError:
                continue
    print(f"{total} windows collected")
