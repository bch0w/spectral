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
    """ Python adaptation of wf_fft.m by Michael West
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

i = 0
path = "/home/bchow/Work/data/egfs/NALASKA_EGF/hyp"

for src in glob(os.path.join(path, "*")):
    src_name = os.path.basename(src)
    for kernel in glob(os.path.join(src, "*")):
        kernel_name = os.path.basename(kernel)
        i = 0
        amp_arr = None
        for fid in glob(os.path.join(kernel, "*")):
            st = read(fid)
            tr = st[0]  # assuming only one trace per stream

            start = tr.stats.starttime
            end = start + (10 * 60)  # 10 minutes max given longest src-rcv distance
            tr.trim(start, end)

            nyquist = tr.stats.sampling_rate / 2
            
            amp, phase, freq = wf_fft(tr.data, nyquist)
            if amp_arr is None:
                amp_arr = amp
            else:
                amp_arr += amp

            i += 1

        plt.plot(freq, amp_arr / i, c="k", lw=1.5)

        for period in [5, 8, 10, 20, 50]:
            plt.axvline(1/period, c="r", ls="--")
            plt.text(s=f" {period}s", x=1/period, y=0.5, c="r")

        plt.xlim([0, 0.3])

        plt.xlabel("Freq [Hz]")
        plt.ylabel("Amplitude")
        plt.title(f"{src_name} {kernel_name} {i} Stations")
        plt.savefig(f"{src_name}_{kernel_name}.png")
        plt.close("all")


            

