"""
Plot the harmonic wave equation to understand wave phenomenon
"""
import matplotlib.pyplot as plt
from numpy import pi, cos, sin, arange, radians, linspace


# User set 
A = 1           # amplitude
f = 1           # frequency (Hz)
l = 1           # wavelength (m)
fix = "t"       # which parameter to keep fixed 't'ime or 'x' (space)

# Parameters that follow
w = 2 * pi * f  # angular frequency (Hz)
k = 2 * pi / l  # wavenumber (m)
T = 1 / f       # period (s)

# Defines the harmonic wave equation
def harmonic_wave_eq(x, t):
    u = A * cos(w * t - k * x)
    return u

# For a given fixed point x0, show what amplitudes look like over time
if fix == "x":
    t = linspace(0, pi, 1000)
    for n in arange(0, 4, 1):
        x = n * pi
        y = harmonic_wave_eq(x, t)
        plt.plot(t, y, label=f"x={n}$\pi$")
        plt.xlabel("time [s]")
        plt.title(f"$\omega$=2$\pi${f:.2f}; k={k:.2f}; T={T}")
    # Show period
    plt.plot([0, T], [A, A], "k--")
    plt.legend()
elif fix == "t":
    x = linspace(0, pi, 1000)
    for t in arange(0, 5, 1):
        y = harmonic_wave_eq(x, t)
        plt.plot(x, y, label=f"t={t}")
        plt.xlabel("x [m]")
        plt.title(f"$\omega$=2$\pi${f:.2f}; k={k:.2f}; $\lambda$={l}")
    # Show period
    plt.plot([0, l], [A, A], "k--")
    plt.legend()
plt.show()



