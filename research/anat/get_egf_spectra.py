"""
Calculating spectra of EGF data from Liu et al. (2022) to determine the optimal 
period bands to choose
"""
import os
import matplotlib.pyplot as plt
import numpy as np

from glob import glob
from obspy import read


def wf_fft(wf, fNyq):
    """ 
    Python adaptation of wf_fft.m by Michael West
    Necessary for GEOS626 work

    INPUT:
        wf - Numpy array of the data points in your trace
        fNyq - the Nyquist frequency

    OUTPUT:
        fft_amp - Numpy array of spectral amplitudes
        fft_phase - Numpy array of phases
        fft_freq - Numpy array of frequencies"""

    NFFT = int(2 ** (np.ceil(np.log(len(wf)) / np.log(2))))  # Next highest power of 2
    FFTX = np.fft.fft(wf, n=NFFT)  # Take fft, padding with zeros.
    NumUniquePts = int(np.ceil((NFFT + 1) / 2))
    FFTX = FFTX[0:NumUniquePts]  # throw out neg frequencies
    MX = abs(FFTX)  # Take magnitude of X
    MX = MX * 2  # Multiply by 2
    fft_amp = MX / len(wf)

    fft_phase = np.angle(FFTX)  # Take magnitude of X

    f = (np.arange(NumUniquePts)) * 2 / NFFT
    fft_freq = f * fNyq

    return fft_amp, fft_phase, fft_freq


if __name__ == "__main__":
    j = 0
    # File structure path/to/<src_station>/ZZ/*.SAC 
    path = "./"

    for src in glob(os.path.join(path, "*")):
        src_name = os.path.basename(src)

        for kernel_name in ["ZZ", "TT"]:
            amp_arr = []
            for fid in glob(os.path.join(src, kernel_name, "*")):
                st = read(fid)
                tr = st[0]  # assuming only one trace per stream

                start = tr.stats.starttime
                # 10 minutes max given longest src-rcv distance
                end = start + (10 * 60)  
                tr.trim(start, end)

                nyquist = tr.stats.sampling_rate / 2
                
                amp, phase, freq = wf_fft(tr.data, nyquist)
                amp_arr.append(amp)

            # Get mean and 1std from the amp array
            amp_arr = np.array(amp_arr)
            mean_arr = amp_arr.T.mean(axis=1)
            std_arr = amp_arr.T.std(axis=1)

            # Plot the data
            f, ax = plt.subplots(dpi=200)

            for amp in amp_arr:
                plt.plot(freq, amp, c="k", alpha=0.05, zorder=3)

            plt.plot(freq, mean_arr, c="k", lw=1.5, zorder=8)
            # plt.plot(freq, mean_arr + std_arr, c="k", lw=1, ls="--", alpha=0.75)
            # plt.plot(freq, mean_arr - std_arr, c="k", lw=1, ls="--", alpha=0.75)

            for period in [1, 5, 8, 10, 20, 30, 100]:
                plt.axvline(1/period, c="r", ls="-", alpha=0.5, zorder=5)
                plt.text(s=f" {period}s", x=1/period, y=0.5, c="r", 
                        fontsize=8, zorder=10)

            # plt.yscale("log")
            plt.xlim([0, 0.3])
            plt.ylim([0, mean_arr.max()])

            plt.xlabel("Freq [Hz]")
            plt.ylabel("Amplitude")
            plt.title(f"{src_name} {kernel_name} {len(amp_arr)} Stations")
            plt.savefig(f"figures/{src_name}_{kernel_name}.png")

            plt.close("all")

