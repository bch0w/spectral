import os
import glob
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from pyatoa.utils.operations.conversions import ascii_to_mseed


def quick_preprocess(st_in):
    """quickly preprocess a stream object"""
    st_out = st_in.copy()
    st_out.detrend("demean")
    st_out.detrend("linear")
    st_out.taper(max_percentage=0.05)
    st_out.filter("bandpass", freqmin=1/30, freqmax=1/10)
    st_out.taper(max_percentage=0.05)

    return st_out


if __name__ == "__main__":
    path_a = "./nz_north_benchmark"
    path_b = "./nz_north_checker_hann_80km"

    paths = [os.path.basename(_) for _ in [path_a, path_b]]
    colors = ['b', 'g']

    for fid_a in glob.glob(os.path.join(path_a, "*semd")):
        if "NZ" not in os.path.basename(fid_a):
            continue
        fid_b = os.path.join(path_b, os.path.basename(fid_a))
        if os.path.exists(fid_b):
            for i, fid in enumerate([fid_a, fid_b]):
                st = ascii_to_mseed(fid, UTCDateTime('2000-01-01T00:00:00'))
                st = quick_preprocess(st)
                plt.plot(st[0].data, c=colors[i], linewidth=2, label=paths[i])

            plt.title("{s}\n{a} vs {b}\n[10,30]s".format(s=os.path.basename(fid_a),
                                              a=paths[0], b=paths[1])
                      )
            plt.grid(True)
            plt.legend()
            plt.xlabel("Samples")
            plt.ylabel("Amplitude [m]")
            plt.savefig(os.path.basename(fid_a) + ".png")
            # plt.show()
            plt.close()



