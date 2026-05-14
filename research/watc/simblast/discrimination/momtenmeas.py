"""
Moment Tensor Measurements (MomTenMeas)

Get Instaseis synthetics for S-wave amplitudes for various moment tensors along
the lune to determine the emergence of S-waves based on deviation from an 
isotropic source. The expectation here is that there is no S-wave energy for an 
isotropic source but introduction of a DC or CLVD component will start to 
generate S-waves at the source. The question is at what rate and at what pct
DC or CLVD are S-waves significant enough to affect classic body-wave 
discrimination techniques
"""
import os
import sys
import instaseis
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pygmt
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from matplotlib.patches import Rectangle
from obspy import read, read_events, Stream
from obspy.geodetics import kilometers2degrees, gps2dist_azimuth
from obspy.taup import TauPyModel

# Import plotting from the internal spectral Class
import sys
sys.path.insert(0, "/Users/prof/Repos/spectral")
from plotools.prettyplot import PrettyPlot


# For TauP Phase list
P_TRAIN = ["p", "P", "PP", "pP", "Pn", "Pg"]
S_TRAIN = ["s", "S", "SS", "sS", "Sn", "Sg"]  


# Toggles
SKIP = False
SHOW = True


def cmaphex(nvals, cmap="seismic"):
    """Return a list of hex codes for `nvals` of a `cmap`"""
    cmap = plt.get_cmap(cmap, nvals)
    hex_out = []
    for i in range(cmap.N):
        rgba = cmap(i)
        # rgb2hex accepts rgb or rgba
        hex_out.append(mpl.colors.rgb2hex(rgba))
    return hex_out


class MomTenMeas:
    """
    For a given moment tensor, query Instaseis/AxiSEM synthetics and make
    amplitude measurements that can be used for discrimination purposes
    """
    def __init__(self, dist_km, baz, src_depth_km, tmin, tmax, choice="ps", 
                 p_phase_list=None, s_phase_list=None, arrival_choice="taup",
                 components="ZNE", kind="velocity", syngine="iasp91_2s",  
                 taup_model="iasp91", taup_buffer=0, 
                 fig_path="FIGURES", wav_path="SAC"):
        """
        Set up calling structure
        
        Parameters
        ----------
        dist_km : float
            Distance from source to receiver in kilometers
        baz : float
            Back azimuth from source to receiver in degrees
        src_depth_km : float
            Source depth in kilometers
        tmin : float
            Minimum time window in seconds
        tmax : float
            Maximum time window in seconds
        choice : str, optional
            Phase choice for measurement ('ps', 'p', 's'), default is "ps"
        p_phase_list : list, optional
            List of P-wave phase names to consider in TauP
        s_phase_list : list, optional
            List of S-wave phase names to consider in TauP
        arrival_choice : str, optional
            Method for determining phase arrivals 
            - 'taup': Use TauP with with `taup_model`
            - 'group': Select based on group arrival wavespeeds
        components : str, optional
            Seismogram components to extract ('ZNE'), default is "ZNE"
        kind : str, optional
            Waveform type ('velocity', 'displacement', etc.), default is "velocity"
        syngine : str, optional
            Syngine model to use for synthetics, default is "iasp91_2s"
        taup_model : str, optional
            TauP velocity model, default is "iasp91"
        taup_buffer: float, optional
            TauP usually provides single travel times. This option pads the TT
            by a given percentage so that each arrival constitutes a window. 
            Allows for some variation with longer period arrivals.
        fig_path : str, optional
            Directory path to save figures, default is "FIGURES"
        wav_path : str, optional
            Directory path to save waveform data, default is "SAC"
        """
        # Input parameters for processing
        self.dist_km = dist_km
        self.dist_deg = kilometers2degrees(kilometer=self.dist_km)

        self.baz = baz
        self.src_depth_km = src_depth_km
        self.tmin = tmin
        self.tmax = tmax
        self.choice = choice
        self.p_phase_list = p_phase_list or P_TRAIN
        self.s_phase_list = s_phase_list or S_TRAIN
        self.arrival_choice = arrival_choice
        assert(self.arrival_choice in ["taup", "group", "custom"])
        self.components = components
        self.kind = kind

        # Set up the Instaseis Database
        self.taup_model = taup_model
        self.syngine = syngine
        self.db = instaseis.open_db(f"syngine://{self.syngine}")
        self.taup_model = taup_model
        self.taup_buffer = taup_buffer
        
        # For storing output files
        self.fig_path = fig_path
        self.wav_path = wav_path
        for path_ in [self.fig_path, self.wav_path]:
            if not os.path.exists(path_):
                os.makedirs(path_)

        # Empty vals for later
        self.idx = None
        self.rcv = None
        self.src = None
        self.st = None
        self.pwin_s = None
        self.swin_s = None
        self.meas = None
        self.save_tag = None

    def get_src(self, latitude=0, longitude=0):
        """
        Turn one of the event objects (based on idx) into an InstaSeis src 
        object for synthetic lookup
        """
        mt = self.event.focal_mechanisms[0].moment_tensor.tensor
        return instaseis.Source(
            latitude=latitude, longitude=longitude, 
            depth_in_m=self.src_depth_km * 1E3,
            m_rr=mt["m_rr"], m_tt=mt["m_tt"], m_pp=mt["m_pp"], 
            m_rt=mt["m_rt"], m_rp=mt["m_rp"], m_tp=mt["m_tp"]
            )
    
    def get_rcv(self, network="XX", station=None):
        """
        Define a receiver lat/lon location based on backazimuth and distance and
        return the Instaseis receiver object
        """
        if station is None:
            station = f"{self.idx:0>4}"
        # Very rudimentary calculation in Cartesian coordinates, shooting 
        # down a backazimuth for a given distance
        lon = self.dist_deg * np.cos(np.deg2rad(self.baz))
        lat = self.dist_deg * np.sin(np.deg2rad(self.baz))

        lon += self.src.latitude  # offset by event lat
        lon += self.src.longitude  

        rcv = instaseis.Receiver(latitude=lat, longitude=lon,
                                 network=network, station=station)      
        return rcv

    def get_synthetics(self):
        """
        Query Instaseis Database for synthetic waveforms, preprocess and attach 
        SAC headers. Process and save waveforms, too. We filter the waveforms 
        before saving so that when we plot with PrettyPlot we get the exact
        same waveforms, instead of having to repeat and possibly change the
        processing steps taken.
        """
        st = self.db.get_seismograms(source=self.src, receiver=self.rcv, 
                                     components=self.components,  
                                     kind=self.kind, dt=0.1)

        # Calcualte dist, az and baz from the now set locations
        dist_m, az, baz = gps2dist_azimuth(lat1=self.src.latitude, 
                                         lon1=self.src.longitude,
                                         lat2=self.rcv.latitude,
                                         lon2=self.rcv.longitude)
        dist_km = dist_m * 1E-3
        dist_deg = kilometers2degrees(dist_km)
        
        # Append SAC Header, copied from PySEP
        sac_header = {
            "iztype": 9,  # Ref time equivalence, IB (9): Begin time
            "b": st[0].stats.starttime,
            "e": st[0].stats.endtime,
            "evla": self.src.latitude,
            "evlo": self.src.longitude,
            "evdp": self.src.depth_in_m * 1E-3,
            "stla": self.rcv.latitude,
            "stlo": self.rcv.longitude,
            "stel": 0,
            "kevnm": str(self.idx),
            "nzyear": st[0].stats.starttime.year,
            "nzjday": st[0].stats.starttime.julday,
            "nzhour": st[0].stats.starttime.hour,
            "nzmin": st[0].stats.starttime.minute,
            "nzsec": st[0].stats.starttime.second,
            "nzmsec": int(f"{st[0].stats.starttime.microsecond:0>6}"[:3]),  # micros->ms see #152
            "dist": dist_km,
            "az": az,  # degrees
            "baz": baz,  # degrees
            "gcarc": dist_deg,  # degrees
            "lpspol": 0,  # 1 if left-hand polarity (usually no in passive seis)
            "lcalda": 1,  # 1 if DIST, AZ, BAZ, GCARC to be calc'd from metadata
        }

        for tr in st:
            tr.stats.sac = sac_header
        
        # Write out SAC File for the RAW waveform
        for tr in st:
            fid = f"{self.wav_path}/{self.save_tag}_{tr.stats.component}.SAC"
            tr.write(fid, format="SAC")

        return st
    
    def load_synthetics(self):
        """
        Load synthetics from disk to avoid re-querying Instaseis each time.
        Expectation is that synthetics were written by `get_synthetics` as SAC
        files and will have already been processed
        """
        st = Stream()
        for comp in self.components:
            fid = f"{self.wav_path}/{self.save_tag}_{comp}.SAC"
            if os.path.exists(fid):
                st += read(fid)
            else:
                raise FileNotFoundError
        return st

    def get_taup_arrivals(self, p_phase_list=None, s_phase_list=None):
        """
        Get arrival time windows from TauP for a given `TAUP_MODEL`. Returns 
        expected P and S-wave arrival windows in units of `samples`
        """
        model = TauPyModel(model=self.taup_model)
        p_phase_list = p_phase_list or self.p_phase_list
        s_phase_list = s_phase_list or self.s_phase_list

        # Get P-phase windows
        p_arrivals = model.get_travel_times(
            source_depth_in_km=self.src_depth_km, 
            distance_in_degree=self.dist_deg,
            phase_list=p_phase_list
            )
        if not p_arrivals:
            raise Exception(f"No P-arrivals found for {self.idx}")

        # Take only time information, add a buffer around direct pick 
        p_arrivals = [_.time for _ in p_arrivals]
        p_window = [min(p_arrivals) * (1-self.taup_buffer) - self.tmax, 
                    max(p_arrivals) * (1+self.taup_buffer) + self.tmax]

        # Get S-phase windows
        s_arrivals = model.get_travel_times(
            source_depth_in_km=self.src_depth_km, 
            distance_in_degree=self.dist_deg,
            phase_list=s_phase_list
            )
        if not s_arrivals:
            raise Exception(f"No P-arrivals found for {self.idx}")
        
        # Take only time information, add a buffer around direct pick and 
        # enlarge the buffer by the largest period
        s_arrivals = [_.time for _ in s_arrivals]
        s_window = [min(s_arrivals) * (1-self.taup_buffer) - self.tmax , 
                    max(s_arrivals) * (1+self.taup_buffer) + self.tmax]

        return p_window, s_window
    
    def get_group_vel_arrivals(self):
        """
        Return phase arrivals based on regional India phase velocities provided
        by Rodgers et al. 2001
        """
        group_velocities = {
            "Pn": [7.7, 8.25],  # km/s
            "Pg": [5.5, 6.5],
            "Sn": [4.0, 4.6],
            "Sg": [3.0, 3.6]  # Lg, values from Baker et al. 2012
            # "p_train": ,
            # "s_train": [6],
        }
        
        # Figure out arrival time based on straight line distance
        arrivals = {}
        for key, vals in group_velocities.items():
            # [min_arrival, max_arrival]
            arrivals[f"{key}"] = [self.dist_km / val for val in vals][::-1]  # s

        p_window = [np.inf, -np.inf]
        s_window = [np.inf, -np.inf]
        for phase, arvs in arrivals.items():
            # Go through and figure out min and max arrival times based on phase
            if phase in self.p_phase_list:
                if min(arvs) < p_window[0]:
                    p_window[0] = min(arvs)
                if max(arvs) > p_window[1]:
                    p_window[1] = max(arvs)
            if phase in self.s_phase_list:
                if min(arvs) < s_window[0]:
                    s_window[0] = min(arvs)
                if max(arvs) > s_window[1]:
                    s_window[1] = max(arvs)

        return p_window, s_window
    
    def get_arrivals_custom(self):
        """
        Custom windows based on picking from waveform plots for 2-4Hz waveforms
        """
        dict_out = {
            150: [[20, 40], [41.75, 55]],
            250: [[35, 44.5], [67, 77]],
            500: [[], []],
            1000: [[120, 180], [233, 297]],
        }
        p_win, s_win = dict_out[self.dist_km]

        return p_win, s_win
    
    def make_measurement(self):
        """
        Take waveforms and calculate a certain measurement, return the 
        measurement value

        Parameters
        - choice (str): 'ps'=P/S amplitude ratio, 's': maximum S amplitude,
            'p': maximum P amplitude
        """
        # Convert from units time to sampling rate for plotting and picking
        samprate = self.st[0].stats.sampling_rate

        # Figure out max amplitude and corresponding index within P window
        p_start, p_end = [int(_ * samprate) for _ in self.pwin_s]  # unit: samp
        s_start, s_end = [int(_ * samprate) for _ in self.swin_s]       
        maxdict = {}
        p_max_avg, s_max_avg = 0, 0
        for tr in self.st:
            # Segment the data in the given window
            p_data_win  = tr.data[p_start:p_end]
            s_data_win  = tr.data[s_start:s_end]

            # Determine global absolute maximum value
            p_max_idx = np.argmax(np.abs(p_data_win)) + p_start
            s_max_idx = np.argmax(np.abs(s_data_win)) + s_start
            
            # Get corresponding max value
            p_max_val = tr.data[p_max_idx]
            s_max_val = tr.data[s_max_idx]    

            # Create an output dictionary for reference later
            maxdict[tr.stats.component] = {"p_max_idx": p_max_idx,
                                           "p_max_val": p_max_val,
                                           "s_max_idx": s_max_idx,
                                           "s_max_val": s_max_val
                                           }
            
            # To calculate the average over all components
            p_max_avg += np.abs(p_max_val)
            s_max_avg += np.abs(s_max_val)

        # Determine the measurement value
        p_max_avg /= len(self.st)
        s_max_avg /= len(self.st)
        ps_ratio = p_max_avg / s_max_avg

        if self.choice.startswith("ps"):
            measurement = ps_ratio
        elif self.choice.startswith("s"):
            measurement = s_max_avg
        elif self.choice.startswith("p"):
            measurement = p_max_avg
        else:
            raise NotImplementedError
        
        # Allow for log10
        if self.choice.endswith("log"):
            raise NotImplementedError("negative values won't plot pygmt")
            measurement = np.log10(measurement)
                
        return maxdict, measurement
        
    def plot_waveforms(self, save="waveform.png", show=False):
        """
        Plot single stream of waveforms with measurement, if available
        """
        # Only plot waveforms for specific moment tensors
        f, axs = plt.subplots(len(self.st), dpi=200, sharex=True)
        # Edge case when we only hve 1 component, still allows for looping
        if len(self.st) == 1:
            axs = [axs]
        middle_idx = len(self.st) // 2

        # Window bounds for plotting rectangles
        pwin_start_s, pwin_end_s = self.pwin_s
        swin_start_s, swin_end_s = self.swin_s

        for i, (tr, ax) in enumerate(zip(self.st, axs)):
            ax.plot(tr.times(), tr.data, label=tr.stats.component, lw=.8, 
                    zorder=6, c="gray")

            maxvals = self.maxdict[tr.stats.component]

            # Plot peak amplitudes as markers
            ax.scatter(maxvals["p_max_idx"] / tr.stats.sampling_rate,
                       maxvals["p_max_val"], c="C0", s=8, edgecolor="k",
                       zorder=7)
            ax.scatter(maxvals["s_max_idx"]  / tr.stats.sampling_rate, 
                       maxvals["s_max_val"], c="C1", s=8, edgecolor="k",
                       zorder=7)

            # Used for scaling
            maxval = np.amax(tr.data)
            
            if i == middle_idx:
                p_label = "P window"
                s_label = "S window"
            else:
                p_label, s_label = None, None

            # Plot windows as transparent rectangels
            p_window_rc = Rectangle(xy=(pwin_start_s, -1.25 * maxval),  
                                    width=pwin_end_s - pwin_start_s,
                                    height=maxval * 2.5, 
                                    alpha=0.25,  color="C0", zorder=2, 
                                    label=p_label)
            s_window_rc = Rectangle(xy=(swin_start_s, -1.25 * maxval),  
                                    width=swin_end_s - swin_start_s,
                                    height=maxval * 2.5, 
                                    alpha=0.25,  color="C1", zorder=2, 
                                    label=s_label)
            ax.add_patch(p_window_rc)
            ax.add_patch(s_window_rc)
            ax.set_ylim([-1.25 * maxvals["s_max_val"], 
                         1.25 * maxvals["s_max_val"]])
            ax.legend(loc="upper right", fontsize=8)

                
        # Plot finalizations
        axs[-1].set_xlabel("Time [s]")
        axs[middle_idx].set_ylabel("Velocity [m/s]")

        title = (f"{self.choice}={self.meas:.2f}; "
                 f"d={self.dist_km:.2f}km baz={self.baz:.2f};\n "
                 f"MT #{self.idx}; {self.syngine}; ")
        if self.tmin and self.tmax:
            title += f"T=[{self.tmin}, {self.tmax}]s"

        axs[0].set_title(title)

        plt.xlim([pwin_start_s * .9, swin_end_s * 1.15])  # cut off long tail
        plt.tight_layout()
        plt.subplots_adjust(hspace=0, wspace=0)
        plt.savefig(save)
        if show:
            plt.show()
        plt.close(f)

    def setup(self, path_cmtsolution):
        """
        Initialize all internal attributes used for making measurements
        Split up from `run` so it can be run separately without going through 
        processing step
        """
        self.event = read_events(path_cmtsolution)[0]

        # e.g.,lune_ipts4_iref5_001
        self.idx = int(self.event.resource_id.id.split("/")[2].split("_")[-1])  
        self.save_tag = f"n{self.idx:0>2}_d{self.dist_km}_b{self.baz}"

        # Get Instaseis synthetics
        self.src = self.get_src()
        self.rcv = self.get_rcv()

        # Get TauP arrivals in units of samples
        if self.arrival_choice == "taup":
            self.pwin_s, self.swin_s = self.get_taup_arrivals()  
        elif self.arrival_choice == "group":
            self.pwin_s, self.swin_s = self.get_group_vel_arrivals()
        elif self.arrival_choice == "custom":
            self.pwin_s, self.swin_s = self.get_arrivals_custom()

    def run(self, path_cmtsolution):
        """
        Main Processing Workflow, for a given CMTSOLUTION and initial state
        setup, get synthetics, taup arrivals, pick peak amplitudes and plot, 
        return information needed for beachball plots.
        """  
        self.setup(path_cmtsolution)

        # Make amplitude measurement
        try:
            self.st = self.load_synthetics()
        except FileNotFoundError:
            self.st = self.get_synthetics()

        if self.tmin and self.tmax:
            self.st.filter("bandpass", freqmin=1/self.tmax, freqmax=1/self.tmin, 
                           zerophase=False)


        self.maxdict, self.meas = self.make_measurement()
        
        self.plot_waveforms(save=f"{self.fig_path}/WAV/wav_{self.save_tag}.png",
                            show=False)
        
        return self.idx, self.meas, self.src.tensor


def mtmrun(dist_km, baz, src_depth_km=1, parallel=True, **kwargs):
    """
    Run `MomTenMeas` for multiple CMTSOLUTIONS, collect the measurement 
    information and then plot
    """
    # Loop through all CMTSOLUTIONS and create waveforms
    path_cmtsolutions = Path("CMTSOLUTION/")
    cmtsolutions = path_cmtsolutions.glob("CMTSOLUTION_*")
    assert(cmtsolutions), f"No source files found for {path_cmtsolutions}"

    # Set up class instance to be used for all runs
    mtm = MomTenMeas(dist_km=dist_km, baz=baz, src_depth_km=src_depth_km,
                     **kwargs)

    idxs, max_amps, tensors = [], [], []
    if parallel:
        with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
            results = executor.map(mtm.run, cmtsolutions)
        for result in results:   
            idx, meas, tensor = result     
            idxs.append(idx)
            max_amps.append(meas)
            tensors.append(tensor)
    else:
        for cmtsolution in cmtsolutions:
            idx, meas, tensor = mtm.run(cmtsolution)
            idxs.append(idx)
            max_amps.append(meas)
            tensors.append(tensor)

    # Kludgy: Run setup on 1 event to run the TauP arrival getter so we have
    # access to the window information which is the same for all sources.
    # Probably a better way to do it but this works.
    mtm.setup(list(path_cmtsolutions.glob("CMTSOLUTION_*"))[0])

    return idxs, max_amps, tensors, mtm.pwin_s, mtm.swin_s


def plot_beachballs(x, y, t, title=None, save=False):
    """
    Plot beachballs based on their assigned index `x` and max amplitude `y` with
    the given moment tensor components `t`
    """
    # Used for coloring beachballs, hex color codes for each index can be used
    # to match waveform plots etc.
    hexvals = cmaphex(nvals=len(t), cmap="viridis")

    # Plot PyGMT Beachball diagrams showing variation of MT with PS ratio
    with pygmt.config(FONT="7.5p"):
        region = [0, max(x) + 1,  min(y) * -1,  max(y) * 1.25]

        projection = "X10c/4c"
        
        # Y label with tick every 
        frame = ["af", "+gwhite", 
                "ya5f1g5+lP/S Amplitude Ratio [ZNE]",
                f"xa1f1g1+l< DC{' '*36}ISO{' '*30}CLVD >",
                ]
        fig = pygmt.Figure()
        fig.basemap(region=region, projection=projection, frame=frame)

        # Title as separate figure becuase I cannot for the life of me figure
        # out how to plot the title with the frame
        fig.text(x=0.2, y=region[-1] * 0.95, text=title, justify="TL")

        for x_, y_, t_ in zip(x, y, t):
            tensor_dict = {
                "mrr": t_[0], "mtt": t_[1], "mff": t_[2],
                "mrt": t_[3], "mrf": t_[4], "mtf": t_[5], 
                "exponent": 1}
            fig.meca(spec=tensor_dict, scale="1.25c", longitude=x_, 
                     latitude=y_, depth=0, compression_fill=hexvals[x_-1], 
                     extension_fill="cornsilk",
                     pen="0.5p,black,solid",)
        
    fig.savefig(save, dpi=500)


def main(dist_km=150, baz=45, src_depth_km=1, tmin=2, tmax=4, components="Z",
         p_phase_list=P_TRAIN, s_phase_list=S_TRAIN, arrival_choice="taup",
         parallel=True, syngine="iasp91_2s", taup_model="iasp91", 
         taup_buffer=0.0, fig_path="FIGURES", wav_path="SAC"):
    """Run and plot"""

    # Used for RS and BB plots
    title = (f"{syngine}, dist={dist_km}km, baz={int(baz%360)}; "
             f"T=[{tmin}, {tmax}]; depth={src_depth_km}km")

    # Main processing workflow
    if not SKIP:
        x, y, t, pwin, swin = mtmrun(dist_km, baz, src_depth_km,
                                     tmin=tmin, tmax=tmax, 
                                     components=components, 
                                     p_phase_list=p_phase_list, 
                                     s_phase_list=s_phase_list, 
                                     arrival_choice=arrival_choice,
                                     syngine=syngine, taup_model=taup_model, 
                                     taup_buffer=taup_buffer,
                                     fig_path=fig_path, wav_path=wav_path, 
                                     parallel=parallel)
        
        # Make beachball plots
        save = f"{fig_path}/bb_d{dist_km}_b{baz}.png"
        plot_beachballs(x, y, t, title, save)

    # Plot record sections
    # Custom look for each of the distances
    customization = {
        150: {"wf_scale": 10, "xlim": [20, 70]}, #, "ylim": [-.25E-3, 1.6E-3]},
        250: {"wf_scale": 10, "xlim": [20, 100]},
        500: {"wf_scale": 30, "xlim": [60, 180]},
        750: {"wf_scale": 40, "xlim": [75, 300]},
        1000: {"wf_scale": 50, "xlim": [100, 370]},
    }
    if SKIP:
        tmarks = None
    else:
        tmarks = pwin + swin

    sac_files = []
    for component in components:
        sac_files += Path(wav_path).glob(f"*_d{dist_km}_b{baz}_{component}.SAC")
    save = f"{fig_path}/rs_d{dist_km}_b{baz}.png"
    kwargs = customization[dist_km]
    pp = PrettyPlot(fids=sac_files, wf_type="recsec", wf_recsec_spacing=5, 
                    fmin=1/tmax, fmax=1/tmin, colors=["viridis"], linewidth=1,
                    ylabel=f"Velocity x {kwargs['wf_scale']} [m/s]",
                    # tp_phases=["ttall"],
                    # tp_phases=P_TRAIN + S_TRAIN,
                    tp_model=taup_model, tp_dist_km=dist_km, 
                    tp_depth=src_depth_km, tp_start=0,
                    tmarks=tmarks, tmarks_c=["C0", "C0", "C1", "C1"], 
                    title=title,  save=save, show=SHOW, legend=True, dpi=200, 
                    transparent=False,
                    **kwargs)
    pp.main()


if __name__ == "__main__":
    # syngine="ak135f_1s"
    main(dist_km=150, baz=45, tmin=1, tmax=4, 
         syngine="ak135f_1s", 
         arrival_choice="custom", src_depth_km=20)

    # for dist_km in [150, 250, 500, 750, 1000]:
    #     if dist_km < 500:
    #         choice = "custom"
    #         tmin = None
    #         tmax = None
    #     else:
    #         choice = "taup"
    #         tmin = 2
    #         tmax = 4
    #     for baz in [0, 45, 89]:
    #         main(dist_km, baz, choice, tmin, tmax)
   
    
