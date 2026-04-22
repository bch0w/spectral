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
import matplotlib.pyplot as plt
import numpy as np
import pygmt
from pathlib import Path
from matplotlib.patches import Rectangle
from obspy import read_events
from obspy.geodetics import kilometers2degrees, gps2dist_azimuth
from obspy.taup import TauPyModel

FIGURES="./FIGURES"


class MomTenMeas:
    """
    For a given moment tensor, query Instaseis/AxiSEM synthetics and make
    amplitude measurements that can be used for discrimination purposes
    """
    def __init__(self, path_cmtsolution="./CMTSOLUTION", syngine="iasp91_2s", 
                 taup_model="iasp91", figure_path="FIGURES"):
        """Set up calling structure"""
        # Set up the Instaseis Database
        self.taup_model = taup_model
        self.syngine = syngine
        self.db = instaseis.open_db(f"syngine://{self.syngine}")

        # Find and check source file
        assert(os.path.exists(path_cmtsolution))
        self.event = read_events(path_cmtsolution)[0]
        # lune_ipts4_iref5_001
        self.tag = int(self.event.resource_id.id.split("/")[2].split("_")[-1])  
        self.taup_model = taup_model
        
        self.figure_path = figure_path
        if not os.path.exists(self.figure_path):
            os.makedirs(figure_path)

        # Empty vals for later
        self.rcv = None
        self.src = None
        self.st = None
        self.pwin = None
        self.swin = None
        self.pmax = None
        self.smax = None
        self.meas = None

    def get_src(self, latitude=0, longitude=0, depth_in_m=1E3):
        """
        Turn one of the event objects (based on idx) into an InstaSeis src 
        object for synthetic lookup
        """
        mt = self.event.focal_mechanisms[0].moment_tensor.tensor
        return instaseis.Source(
            latitude=latitude, longitude=longitude, depth_in_m=depth_in_m,
            m_rr=mt["m_rr"], m_tt=mt["m_tt"], m_pp=mt["m_pp"], 
            m_rt=mt["m_rt"], m_rp=mt["m_rp"], m_tp=mt["m_tp"]
            )
    
    def get_rcv(self, baz=0, dist_km=1E3, network="XX", station="S000"):
        """
        Define a receiver lat/lon location based on backazimuth and distance and
        return the Instaseis receiver object
        """
        # Very rudimentary calculation in Cartesian coordinates, shooting 
        # down a backazimuth for a given distance
        dist_deg = kilometers2degrees(kilometer=dist_km)
        lon = dist_deg * np.cos(np.deg2rad(baz))
        lat = dist_deg * np.sin(np.deg2rad(baz))

        lon += self.src.latitude  # offset by event lat
        lon += self.src.longitude  

        rcv = instaseis.Receiver(latitude=lat, longitude=lon,
                                 network=network, station=station)      
        return rcv

    def get_synthetics(self, components="ZNE", kind="velocity"):
        """
        Query Instaseis Database for synthetic waveforms, preprocess and attach 
        SAC headers. Option to save waveforms    
        """
        st = self.db.get_seismograms(source=self.src, receiver=self.rcv, 
                                     components=components,  kind=kind,
                                     dt=0.1)
        
        # convert from nm/s to m/s  (not needed?)
        # for tr in st:
        #     tr.data *= 1E-9  

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
            "kevnm": self.tag,
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

    def get_taup_arrivals(self, source_depth_in_km, distance_in_km, 
                          p_phase_list=None, s_phase_list=None, model="iasp91",
                          buffer=0.025):
        """
        Get arrival time windows from TauP for a given `TAUP_MODEL`. Returns 
        expected P and S-wave arrival windows in units of `samples`
        """
        # By default query all the crustal directarrivals
        if not p_phase_list:
            # p_phase_list = ["p", "P", "PP", "pP", "Pn", "Pg"]
            p_phase_list = ["p", "P", "Pn", "Pg"]
        if not s_phase_list:
            # s_phase_list = ["s", "S", "SS", "sS", "Sn", "Sg"]
            s_phase_list = ["s", "S", "Sn", "Sg"]

        dist_deg = kilometers2degrees(distance_in_km)

        model = TauPyModel(model=model)
        p_arrivals = model.get_travel_times(
            source_depth_in_km=source_depth_in_km, distance_in_degree=dist_deg,
            phase_list=p_phase_list
            )
        if not p_arrivals:
            return None, None
        
        p_arrivals = [_.time for _ in p_arrivals]
        p_window = [min(p_arrivals) * (1-buffer), max(p_arrivals) * (1+buffer)]

        s_arrivals = model.get_travel_times(
            source_depth_in_km=source_depth_in_km, distance_in_degree=dist_deg,
            phase_list=s_phase_list
            )
        if not s_arrivals:
            return None, None
        
        s_arrivals = [_.time for _ in s_arrivals]
        s_window = [min(s_arrivals) * (1-buffer), max(s_arrivals) * (1+buffer)]

        # Convert from units of 
        samprate = self.st[0].stats.sampling_rate
        pwin = [int(_ * samprate) for _ in p_window]         
        swin = [int(_ * samprate) for _ in s_window]          

        # Incase the window is only 1 timestamp long, make sure we don't have a
        # window len 0
        if pwin[0] == pwin[1]:
            pwin[1] += 1  
        if swin[0] == swin[1]:
            swin[1] += 1  

        return pwin, swin
    
    def make_measurement(self, tmin, tmax, choice="ps"):
        """
        Take waveforms and calculate a certain measurement, return the 
        measurement value

        Parameters
        - choice (str): 'ps'=P/S amplitude ratio, 's': maximum S amplitude,
            'p': maximum P amplitude
        """
        self.st.filter("bandpass", freqmin=1/tmax, freqmax=1/tmin, 
                       zerophase=True)

        # Figure out max amplitude and corresponding index within P window
        p_start, p_end = self.pwin   # unit: samples
        s_start, s_end = self.swin
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

        if choice == "ps":
            measurement = ps_ratio
        elif choice == "s":
            measurement = s_max_avg
        elif choice == "p":
            measurement = p_max_avg
        else:
            raise NotImplementedError
                
        return maxdict, measurement
        
    def plot_waveforms(self, save="waveform.png", title="", show=False):
        """Plot single stream of waveforms with measurement, if available"""
            # Only plot waveforms for specific moment tensors
        f, axs = plt.subplots(len(self.st), dpi=200, sharex=True)
        middle_idx = len(self.st) // 2

        # Convert window bounds to unit seconds to match waveform x-axis
        pwin_start_s, pwin_end_s = [_ / self.st[0].stats.sampling_rate 
                                    for _ in self.pwin]
        swin_start_s, swin_end_s = [_ / self.st[0].stats.sampling_rate 
                                    for _ in self.swin]

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
        axs[0].set_title(title)
        plt.xlim([pwin_start_s * .9, swin_end_s * 1.15])  # cut off long waveform tail
        plt.tight_layout()
        plt.subplots_adjust(hspace=0, wspace=0)
        plt.savefig(save)
        if show:
            plt.show()
        plt.close(f)
    
    def run(self, dist_km, baz, src_depth_m, tmin, tmax, choice="ps",
            components="ZNE"):
        """Main Processing Workflow"""  
        # Get Instaseis synthetics
        self.src = self.get_src(depth_in_m=src_depth_m)
        self.rcv = self.get_rcv(dist_km=dist_km, baz=baz)
        self.st = self.get_synthetics(components=components)

        # Get TauP arrivals in units of samples
        self.pwin, self.swin = self.get_taup_arrivals(
            model=self.taup_model, source_depth_in_km=src_depth_m * 1E-3, 
            distance_in_km=dist_km 
            )
        
        # Make amplitude measurement
        self.maxdict, self.meas = self.make_measurement(tmin=tmin, tmax=tmax, 
                                                        choice=choice)
        self.plot_waveforms(
            title=f"{choice}={self.meas:.2f}; "
                  f"d={dist_km:.2f}km baz={baz:.2f}; "
                  f"T=[{tmin}, {tmax}]s\n"
                  f"MT #{self.tag}; {self.syngine}; ",
                  save=f"{self.figure_path}/{self.tag:0>3}_{dist_km}_{baz}.png",
                  show=False
                  )


def main(dist_km, baz, src_depth_m=1E3, tmin=2, tmax=5):
    """
    Run `MomTenMeas` for multiple CMTSOLUTIONS, collect the measurement 
    information and then plot
    """
    tensors, lune_idxs, max_amps = [], [], []
    path_cmtsolutions = Path("CMTSOLUTION/")
    for fid in sorted(path_cmtsolutions.glob("CMTSOLUTION_*")):
        mtm = MomTenMeas(path_cmtsolution=fid)
        mtm.run(dist_km=dist_km, baz=baz, src_depth_m=src_depth_m, 
                tmin=tmin, tmax=tmax)
        
        tensors.append(mtm.src.tensor)
        lune_idxs.append(mtm.tag)
        max_amps.append(mtm.meas)
    
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
                     latitude=amp, depth=0, compression_fill="red", 
                     extension_fill="cornsilk",
                     pen="0.5p,black,solid",)
        
    fig.savefig(f"{mtm.figure_path}/smg_{int(dist_km)}_{int(baz)}.png", dpi=500)


if __name__ == "__main__":
    for dist_km in [100, 250, 500, 1000]:
        for baz in [0, 45, 89]:
            print(f"{dist_km} {baz}")
            main(dist_km, baz)

    
