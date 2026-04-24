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
import instaseis
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pygmt
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from matplotlib.patches import Rectangle
from obspy import read_events
from obspy.geodetics import kilometers2degrees, gps2dist_azimuth
from obspy.taup import TauPyModel

FIGURES="./FIGURES"

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
                 components="ZNE", kind="velocity", syngine="iasp91_2s", 
                 taup_model="iasp91", fig_path="FIGURES", wav_path="SAC"):
        """Set up calling structure"""
        # Input parameters for processing
        self.dist_km = dist_km
        self.dist_deg = kilometers2degrees(kilometer=self.dist_km)

        self.baz = baz
        self.src_depth_km = src_depth_km
        self.tmin = tmin
        self.tmax = tmax
        self.choice = choice
        self.components = components
        self.kind = kind

        # Set up the Instaseis Database
        self.taup_model = taup_model
        self.syngine = syngine
        self.db = instaseis.open_db(f"syngine://{self.syngine}")
        self.taup_model = taup_model
        
        # For storing output files
        self.fig_path = fig_path
        self.wav_path = wav_path
        for path_ in [self.fig_path, self.wav_path]:
            if not os.path.exists(path_):
                os.makedirs(path_)

        # Empty vals for later
        self.tag = None
        self.rcv = None
        self.src = None
        self.st = None
        self.pwin = None
        self.swin = None
        self.pmax = None
        self.smax = None
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
            station = f"{self.tag:0>4}"
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
        SAC headers. Option to save waveforms    
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
            "kevnm": str(self.tag),
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

        return st

    def get_taup_arrivals(self, p_phase_list=None, s_phase_list=None,
                          buffer=0.025):
        """
        Get arrival time windows from TauP for a given `TAUP_MODEL`. Returns 
        expected P and S-wave arrival windows in units of `samples`
        """
        # By default query all the crustal directarrivals
        if not p_phase_list:
            p_phase_list = ["p", "P", "Pn", "Pg"]
        if not s_phase_list:
            s_phase_list = ["s", "S", "Sn", "Sg"]

        model = TauPyModel(model=self.taup_model)

        # Get P-phase windows
        p_arrivals = model.get_travel_times(
            source_depth_in_km=self.src_depth_km, 
            distance_in_degree=self.dist_deg,
            phase_list=p_phase_list
            )
        if not p_arrivals:
            raise Exception(f"No P-arrivals found for {self.tag}")

        # Take only time information, add a buffer around direct pick 
        p_arrivals = [_.time for _ in p_arrivals]
        p_window = [min(p_arrivals) * (1-buffer), max(p_arrivals) * (1+buffer)]

        # Get S-phase windows
        s_arrivals = model.get_travel_times(
            source_depth_in_km=self.src_depth_km, 
            distance_in_degree=self.dist_deg,
            phase_list=s_phase_list
            )
        if not s_arrivals:
            raise Exception(f"No P-arrivals found for {self.tag}")
        
        # Take only time information, add a buffer around direct pick 
        s_arrivals = [_.time for _ in s_arrivals]
        s_window = [min(s_arrivals) * (1-buffer), max(s_arrivals) * (1+buffer)]

        return p_window, s_window
    
    def make_measurement(self):
        """
        Take waveforms and calculate a certain measurement, return the 
        measurement value

        Parameters
        - choice (str): 'ps'=P/S amplitude ratio, 's': maximum S amplitude,
            'p': maximum P amplitude
        """
        self.st.filter("bandpass", freqmin=1/self.tmax, freqmax=1/self.tmin, 
                       zerophase=True)

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

        if self.choice == "ps":
            measurement = ps_ratio
        elif self.choice == "s":
            measurement = s_max_avg
        elif self.choice == "p":
            measurement = p_max_avg
        else:
            raise NotImplementedError
                
        return maxdict, measurement
        
    def plot_waveforms(self, save="waveform.png", show=False):
        """
        Plot single stream of waveforms with measurement, if available
        """
        # Only plot waveforms for specific moment tensors
        f, axs = plt.subplots(len(self.st), dpi=200, sharex=True)
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
                 f"d={self.dist_km:.2f}km baz={self.baz:.2f}; "
                 f"T=[{self.tmin}, {self.tmax}]s\n"
                 f"MT #{self.tag}; {self.syngine};")
        axs[0].set_title(title)

        plt.xlim([pwin_start_s * .9, swin_end_s * 1.15])  # cut off long tail
        plt.tight_layout()
        plt.subplots_adjust(hspace=0, wspace=0)
        plt.savefig(save)
        if show:
            plt.show()
        plt.close(f)
    
    def run(self, path_cmtsolution):
        """
        Main Processing Workflow, for a given CMTSOLUTION and initial state
        setup, get synthetics, taup arrivals, pick peak amplitudes and plot, 
        return information needed for beachball plots.
        """  
        self.event = read_events(path_cmtsolution)[0]

        # e.g.,lune_ipts4_iref5_001
        self.tag = int(self.event.resource_id.id.split("/")[2].split("_")[-1])  
        self.save_tag = f"n{self.tag:0>3}_d{self.dist_km}_b{self.baz}"

        # Get Instaseis synthetics
        self.src = self.get_src()
        self.rcv = self.get_rcv()
        self.st = self.get_synthetics()

        # Get TauP arrivals in units of samples
        self.pwin_s, self.swin_s = self.get_taup_arrivals()  
        
        # Make amplitude measurement
        self.maxdict, self.meas = self.make_measurement()
        
        # Write out SAC File of the filtered waveform
        for tr in self.st:
            fid = f"{self.wav_path}/{self.save_tag}_{tr.stats.component}.SAC"
            tr.write(fid, format="SAC")

        self.plot_waveforms(save=f"{self.fig_path}/{self.save_tag}.png",
                            show=False)
        
        # Main process prints out the pplot command to run
        if self.tag == 1:

            print(f"pplot *.SAC --tmarks "
                  f"{self.pwin_s[0]:.1f} {self.pwin_s[1]:.1f} "
                  f"{self.swin_s[0]:.1f} {self.swin_s[1]:.1f} "
                  f"--tmark_c r r b b")
        
        return self.src.tensor, self.tag, self.meas


def main(dist_km, baz, src_depth_km=1, tmin=2, tmax=5, parallel=True):
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
                     tmin=tmin, tmax=tmax)

    tensors, lune_idxs, max_amps = [], [], []
    if parallel:
        with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
            results = executor.map(mtm.run, cmtsolutions)
        for result in results:   
            tensor, tag, meas = result     
            tensors.append(tensor)
            lune_idxs.append(tag)
            max_amps.append(meas)
    else:
        for cmtsolution in cmtsolutions:
            tensor, tag, meas = mtm.run(cmtsolution)
            tensors.append(tensor)
            lune_idxs.append(tag)
            max_amps.append(meas)
    
    # Used for coloring beachballs, hex color codes for each index can be used
    # to match waveform plots etc.
    hexvals = cmaphex(nvals=len(tensors), cmap="RdYlBu_r")

    # Plot PyGMT Beachball diagrams showing variation of MT with PS ratio
    with pygmt.config(FONT="7.5p"):
        region = [0, 
                  max(lune_idxs) + 1, 
                  min(max_amps) * -1, 
                  max(max_amps) * 1.25]

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
        fig.text(
            x=0.2, y=region[-1] * 0.95, 
            text=f"{mtm.taup_model}, dist={dist_km:.2f}km, baz={baz%360:.2f}",
            justify="TL")

        for idx, amp, tensor in zip(lune_idxs, max_amps, tensors):
            tensor_dict = {
                "mrr": tensor[0], "mtt": tensor[1], "mff": tensor[2],
                "mrt": tensor[3], "mrf": tensor[4], "mtf": tensor[5], 
                "exponent": 1}
            fig.meca(spec=tensor_dict, scale="1.25c", longitude=idx, 
                     latitude=amp, depth=0, compression_fill=hexvals[idx-1], 
                     extension_fill="cornsilk",
                     pen="0.5p,black,solid",)
        
    fig.savefig(f"{mtm.fig_path}/smg_{int(dist_km)}_{int(baz)}.png", dpi=500)


if __name__ == "__main__":
    main(100, 45, parallel=True)

    # for dist_km in [100, 250, 500, 1000]:
    #     for baz in [0, 45, 89]:
    #         print(f"{dist_km} {baz}")
    #         main(dist_km, baz)

    
