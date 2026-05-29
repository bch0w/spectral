"""
Try to "locate" a known earthquake with small aperture array methods
"""
from obspy import read, Stream
from obspy import UTCDateTime as UTC
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


data_path = Path("/Users/chow/Work/research/gulkanaseis24/nodes")

r_maxs = []
for station in np.arange(100, 112, 1, dtype=int):
    st = Stream()
    for comp in "ZNE":
        st += read(data_path / f"GS.{station}.DH{comp}.2024.253")

    # Cut over the whole waveform
    st.trim(UTC("2024-09-09T01:40:49"), UTC("2024-09-09T01:42:00"))
    st.filter("highpass", freq=1)

    r, t, z, t_sum = [], [], [], []
    bazs = range(0, 360, 1)
    for baz in bazs:
        st_r = st.copy()
        st_r.rotate("NE->RT", back_azimuth=baz)

        # P-wave arrival
        st_r.trim(UTC("2024-09-09T01:40:52.3"), UTC("2024-09-09T01:40:52.44"))

        # Get radial max value
        tr_r = st_r.select(component="R")[0].data
        absmax_r = max(tr_r.min(), tr_r.max(), key=abs)
        r.append(absmax_r)

        # Get sum of transverse
        t_sum.append(np.sum(st_r.select(component="T")[0].data))

    r = np.array(r)
    t_sum = np.array(t_sum)

    r_max = bazs[np.where(r == max(r))[0][0]]
    r_maxs.append(r_max)
    if False:
        plt.plot(bazs, r/r.max(), "r", label=f"absmax radial ({r_max})", alpha=0.75)

        plt.axvline(r_max, c="r", ls="--", alpha=0.75)
        plt.plot(bazs, t_sum/t_sum.max(), "b", label=f"sum transverse", alpha=0.5)
        # plt.axvline(t_max, c="b", ls="--", alpha=0.5)

        plt.legend()
        plt.xlabel("BAz (deg)")
        plt.ylabel("Max Amp [m/s]")
        plt.title("EQ 11882516 BAz Estimator\nEvent Location BAz = 58.25*")
        plt.grid()
        plt.savefig(f"{station}_11882516.png")
        plt.close()

        # Plot waveform for the max radial
        st_r = st.copy()
        st_r.rotate("NE->RT", back_azimuth=r_max)
        st_r.trim(UTC("2024-09-09T01:40:52.3"), UTC("2024-09-09T01:40:52.44"))
        st_r.plot(show=False, outfile=f"{station}_11882516_wf.png")

# Make a rose plot
plt.hist(r_maxs, bins=np.arange(0, 360, 10))
plt.grid()
plt.xlabel("BAz (deg)")
plt.ylabel("Count")
plt.title("Estimated BAz 11882516")

if False:

    # Trying to get P-wave velocities
    p_picks = {
            "100": UTC("2024-09-09T01:40:52.37586"),
            "101": UTC("2024-09-09T01:40:52.37598"),
            "102": UTC("2024-09-09T01:40:52.37599"),
            "103": UTC("2024-09-09T01:40:52.37599"),
            "104": UTC("2024-09-09T01:40:52.37600"),
            "105": UTC("2024-09-09T01:40:52.37999"),
            "106": UTC("2024-09-09T01:40:52.38402"),
            "107": UTC("2024-09-09T01:40:52.38042"),
            "108": UTC("2024-09-09T01:40:52.37994"),
            "109": UTC("2024-09-09T01:40:52.39598"),
            "110": UTC("2024-09-09T01:40:52.38794"),
            "111": UTC("2024-09-09T01:40:52.39200"),
            }

    s_picks = {
            "100": UTC("2024-09-09T01:41:11.89588"),
            "101": UTC("2024-09-09T01:41:11.89603"),
            "102": UTC("2024-09-09T01:41:11.89587"),
            "103": UTC("2024-09-09T01:41:11.89993"),
            "104": UTC("2024-09-09T01:41:11.90002"),
            "105": UTC("2024-09-09T01:41:11.89596"),
            "106": UTC("2024-09-09T01:41:11.90796"),
            "107": UTC("2024-09-09T01:41:11.90780"),
            "108": UTC("2024-09-09T01:41:11.89977"),
            "109": UTC("2024-09-09T01:41:11.91195"),
            "110": UTC("2024-09-09T01:41:11.91193"),
            "111": UTC("2024-09-09T01:41:11.90791"),
            }

    p_timediffs = {
        "100": 0.0,
        "101": -0.00012,
        "102": -0.00013,
        "103": -0.00013,
        "104": -0.00014,
        "105": -0.00413,
        "106": -0.00816,
        "107": -0.00456,
        "108": -0.00408,
        "109": -0.02012,
        "110": -0.01208,
        "111": -0.01614,
    }

    s_timediffs = {
        "100": 0.000000,
        "101": -0.000150,
        "102": 0.000010,
        "103": -0.004050,
        "104": -0.004140,
        "105": -0.000080,
        "106": -0.012080,
        "107": -0.011920,
        "108": -0.003890,
        "109": -0.016070,
        "110": -0.016050,
        "111": -0.012030,
    }


    locations = {
        "100": (63.265612, -145.415025),
        "101": (63.265362, -145.41407),
        "102": (63.265448, -145.41381),
        "103": (63.265203, -145.415355),
        "104": (63.265062, -145.414323),
        "105": (63.26511, -145.413887),
        "106": (63.26495, -145.415472),
        "107": (63.264683, -145.414657),
        "108": (63.264878, -145.414107),
        "109": (63.26392, -145.41601),
        "110": (63.26466, -145.41557),
        "111": (63.264278, -145.415858),
    }

    dist_m = {
        "100": 0.00,
        "101": 55.46,
        "102": 63.69,
        "103": 48.51,
        "104": 70.72,
        "105": 79.98,
        "106": 77.13,
        "107": 105.19,
        "108": 93.91,
        "109": 194.98,
        "110": 109.59,
        "111": 154.47,
    }
