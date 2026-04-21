"""
Get Instaseis synthetics for S-wave amplitudes for various moment tensors along
the lune to determine the emergence of S-waves based on deviation from an 
isotropic source. The expectation here is that there is no S-wave energy for an 
isotropic source but introduction of a DC or CLVD component will start to 
generate S-waves at the source. The question is at what rate and at what pct
DC or CLVD are S-waves significant enough to affect classic body-wave 
discrimination techniques
"""
import instaseis
import matplotlib.pyplot as plt

from pathlib import Path
from matplotlib.patches import Rectangle
from obspy.geodetics import kilometers2degrees, gps2dist_azimuth
from obspy.imaging.beachball import beach
from obspy.taup import TauPyModel


TAUP_MODEL = "iasp91"
SYNGINE = "iasp91_2s"
LOWPASS_S = 2
HIGHPASS_S = 5
P_PHASE_LIST = ["p", "P", "PP", "pP", "Pn", "Pg"]
S_PHASE_LIST = ["s", "S", "SS", "sS", "Sn", "Sg"]


def get_taup_arrivals(source_depth_in_km, distance_in_km, p_phase_list=None,
                      s_phase_list=None, model=TAUP_MODEL):
    """
    Get arrival time windows from TauP for a given `TAUP_MODEL`
    """
    # By default query all the crustal directarrivals
    if not p_phase_list:
        p_phase_list = P_PHASE_LIST
    if not s_phase_list:
        s_phase_list = S_PHASE_LIST

    dist_deg = kilometers2degrees(distance_in_km)

    model = TauPyModel(model=model)
    p_arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km,
                                        distance_in_degree=dist_deg,
                                        phase_list=p_phase_list)
    if not p_arrivals:
        return None, None
    
    p_arrivals = [_.time for _ in p_arrivals]
    p_window = [min(p_arrivals), max(p_arrivals)]

    s_arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km,
                                        distance_in_degree=dist_deg,
                                        phase_list=s_phase_list)
    if not s_arrivals:
        return None, None
    
    s_arrivals = [_.time for _ in s_arrivals]
    s_window = [min(s_arrivals), max(s_arrivals)]

    return p_window, s_window


def cmtsolution_to_instaseis_source(fid, latitude=0., longitude=0.,
                                    depth_in_m=1000.):
    with open(fid, "r") as f:
        lines = f.readlines()
    vals = [float(line.split(":")[1].strip()) for line in lines[7:]]  
    m_rr, m_tt, m_pp, m_rt, m_rp, m_tp = vals

    src = instaseis.Source(
        latitude=latitude, longitude=longitude, depth_in_m=depth_in_m,
        m_rr=m_rr, m_tt=m_tt, m_pp=m_pp, m_rt=m_rt, m_rp=m_rp, m_tp=m_tp
        )
    return src


def receivers_from_stations_file(fid="STATIONS", select=None):
    """
    Get instaseis Receiver objects for all available stations in file
    """
    with open(fid, "r") as f:
        lines = f.readlines()
    receivers = []
    for line in lines:
        sta, net, lat, lon, *_ = line.split()
        if select and sta not in select:
            continue
        receivers.append(instaseis.Receiver(latitude=lat, longitude=lon,
                                            network=net, station=sta))

    return receivers

    
def get_measurement(db, src, rcv, lune_idx, plot=False, choice="s_max"):
    """
    Get the maximum amplitude of an instaseis synthetic for a given window of 
    time calculated using TauP using the same model

    choice (str): `s_max`= S amplitude, `p_max`=P amplitude, 
        `ps_ratio`=P/S ratio
    """
    st = db.get_seismograms(source=src, receiver=rcv, components="ZNE", 
                            kind="velocity")
    for tr in st:
        tr.data *= 1E-9  # convert from nm/s to m/s
    st.filter("bandpass", freqmin=1/HIGHPASS_S, freqmax=1/LOWPASS_S, 
              zerophase=True)
    
    # Get src-rcv distance
    dist_m, az, baz = gps2dist_azimuth(src.latitude, src.longitude,
                                       rcv.latitude, rcv.longitude)
    dist_km = dist_m * 1E-3

    # Determine S-arrival window from TauPy
    p_window, s_window = get_taup_arrivals(
        source_depth_in_km=src.depth_in_m * 1E-3, 
        distance_in_km=dist_km, p_phase_list=P_PHASE_LIST, 
        s_phase_list=S_PHASE_LIST, 
        model=TAUP_MODEL
        )

    # window start and end based on samples
    sr = st[0].stats.sampling_rate
    p_start_idx, p_end_idx = [int(_ * sr) for _ in p_window]          
    s_start_idx, s_end_idx = [int(_ * sr) for _ in s_window]          

    # Incase the window is only 1 timestamp long
    if p_start_idx == p_end_idx:
        p_end_idx += 1  
    if s_start_idx == s_end_idx:
        s_end_idx += 1  

    # Take the largest amplitude from each, assuming it is the S-wave
    # positive values only
    p_max_vals = [abs(tr.data[p_start_idx:p_end_idx].max()) for tr in st]
    s_max_vals = [abs(tr.data[s_start_idx:s_end_idx].max()) for tr in st]

    s_max_single = max(s_max_vals)  # for scaling the plot of waveforms

    p_max = sum(p_max_vals) / 3
    s_max = sum(s_max_vals) / 3
    ps_ratio = p_max / s_max

    if choice == "pmax":
        return_val = p_max
    elif choice == "smax":
        return_val = s_max
    elif choice == "ps_ratio":
        return_val = ps_ratio

    # Only plot waveforms for specific moment tensors
    if plot and (lune_idx in [1, 8, 15]):
        f, ax = plt.subplots(1, dpi=200)

        for i, tr in enumerate(st):
            # Scale by max S amplitude for better visualization of 
            # P and S arrivals
            _p_max = p_max_vals[i]  
            _s_max = s_max_vals[i]
            plt.plot(tr.times(), tr.data / s_max_single + i, 
                     label=(f"{tr.stats.component}"
                            f"(P_max={_p_max:.2E}m/s; S_max={_s_max:.2E}m/s)"),
                     lw=1, zorder=6)
            
            # Plot windows for P and S arrival prediction windows
            ymin, ymax = plt.gca().get_ylim()

        labelled = False
        for window, label in zip([p_window, s_window], 
                                 ["P window", "S window"]):
            if not labelled:
                label = label
                labelled = True
            rc = Rectangle(
                xy=(window[0],  ymin),  width=window[1]-window[0],
                height=ymax-ymin, alpha=0.25, 
                color=f"C{i+3}", zorder=2, label=label
                )
            i += 1  # just to get a different color for the S window 
            plt.gca().add_patch(rc)
                
        plt.xlim([0, s_window[1] * 1.25])
        plt.ylim([-2, 4])
        plt.xlabel("Time [s]")
        plt.ylabel("Normalized Amplitude")
        plt.title(f"LUNE_{lune_idx} {db.model}; "
                f"dist={dist_km:.2f}km; baz={baz:2f}",)
        plt.legend(loc="upper right", fontsize=8)
        plt.tight_layout()
        plt.savefig(plot)
        plt.close(f)

    return dist_km, baz, return_val


def main(dist):
    """
    Open Syngine database, calculate maximum amplitude for a given station 
    epicentral distance for all available moment tensors. Plot on a scatter 
    plot with indices corresponding to the lune indices from Carl's plot
    """
    instaseis_db = f"syngine://{SYNGINE}"
    rcv_lat = dist/111.111
    rcv_lon = 0
    choice = "ps_ratio"
    threshold = 0.8

    cmtsolutions = Path("CMTSOLUTION")
    cmtsolutions = sorted(list(cmtsolutions.glob("CMTSOLUTION_*")))

    assert(bool(len(cmtsolutions))), "no source files found" 

    db = instaseis.open_db(instaseis_db)
    rcv = instaseis.Receiver(latitude=rcv_lat, longitude=rcv_lon,
                             network="NN", station="S000")

    f, ax = plt.subplots(1, dpi=200)
    tensors, max_amps, lune_idxs = [], [], []
    for i, cmtsolution in enumerate(cmtsolutions):
        src = cmtsolution_to_instaseis_source(cmtsolution)
        lune_idx = int(cmtsolution.name.split("_")[-1])

        dist_km, baz, max_amp = get_measurement(
            db, src, rcv, lune_idx=lune_idx, choice=choice,
            plot=f"FIGURES/WAV/wav_{dist:.0f}_{lune_idx:0>2}.png"
            )

        # Collect for later scatterplot to connect all MTs
        tensors.append(src.tensor)
        max_amps.append(max_amp)
        lune_idxs.append(lune_idx)

    # Used for scaling positions on figure
    scale = max(max_amps) - min(max_amps)
    pct = 0.0075 * scale
    width = scale / 10

    # Plot all with a scatterplot to set up figure
    plt.plot(lune_idxs, max_amps, "ko-", zorder=10)

    plt.axhline(threshold, c="k", ls="--", zorder=5, alpha=0.5)

    # Plot all beachballs
    for tensor, lune_idx, max_amp in zip(tensors, lune_idxs, max_amps):
        # plot the beachball at the first point
        # Mrr=Mzz, Mtt=Mxx, Mpp=Myy, Mrt=Mxz, Mrp=-Myz, Mtp=-Mxy
        # 1,2,3 = Up,South,East which equals r,theta,phi
        # NM x 6 (M11, M22, M33, M12, M13, M23
        # Mrr, Mtt, Mpp, Mrt, Mrp, Mtp
        bb = beach(tensor, xy=(lune_idx, max_amp), 
                   width=(1, width), facecolor="r", edgecolor="k", 
                   zorder=11, linewidth=1)
        ax.add_collection(bb)
        plt.text(lune_idx, max_amp-pct*10, lune_idx, c="k", zorder=12, 
                 ha="center")

    # Set xticks lines for every lune index [1, n]
    if False:
        xvals = [_ + 1 for _ in list(range(0, i+1))]
        plt.xticks(xvals)
    else:
        # No xticks, rely on text labels next to each beachball
        plt.tick_params(axis="x", which="both", bottom=False, top=False, 
                        labelbottom=False)

    # These indices correspond to the lune MTs that denote a transition from one
    # set to another to help clarify where on the lune we are
    labels = ["DC", "ISO", "CLVD", "DC"]
    idxs = [1, 8, 15, 20.5]
    for lab, val in zip(labels, idxs):
        plt.axvline(val, c="r", ls="-", zorder=8)
        plt.text(val, max(max_amps) + 5 * pct, lab, rotation=90, c="r", 
                 zorder=9, ha="left")

    plt.ylim([min(max_amps) - 15 * pct, max(max_amps) + 20 * pct])  # add some buffer for text
    plt.xlabel("Lune Index")
    if choice == "p_max":
        plt.ylabel("Max S-wave Amplitude [ZNE]")
    elif choice == "s_max":
        plt.ylabel("Max S-wave Amplitude [ZNE]")
    elif choice == "ps_ratio":
        plt.ylabel("P/S Amplitude Ratio [ZNE]")
    plt.title(f"{TAUP_MODEL.upper()}; dist={dist_km:.2f}km; BAz={baz%360:.2f}")
    plt.tight_layout()
    plt.savefig(f"FIGURES/smg_{dist:.0f}.png")
    plt.close(f)


if __name__ == "__main__":
    for dist in [50, 100, 250, 500, 1000, 1500, 2500, 5000]:
        main(dist=dist)
    
