"""
Plotting Rayleigh wave sensitivity for homogeneous isotropic Poisson's solid

TO DO: Plot arrows at depth layers representing u_z and u_x to show retrograde
to prograde transition
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import e, pi, sin, cos


T = 1000
vs = 3.5  # S-wave from IASP91 [km/s]
vr = vs * 0.92

lambda_x = vr * T  # km

# lambda_x = 1
k_x = 2 * pi / lambda_x

z = np.linspace(0, lambda_x * 2.5, 1000)
u_x =  e**(-0.85 * k_x * z) - 0.58 * e **(-0.39 * k_x * z)
u_z =  -0.85 * e**(-0.85 * k_x * z) + 1.47 * e **(-0.39 * k_x * z)

f, ax = plt.subplots(figsize=(4.5,8))
plt.axvline(0, ls="--", lw=1, c="k")
plt.plot(u_x, z, c="C0", label="ux")
plt.plot(u_z, z, c="C1", label="uz")
plt.ylim([0, z.max()])
ax.invert_yaxis()

plt.legend()
plt.grid()
plt.xlabel("displacement")
if lambda_x == 1:
    plt.ylabel("depth/$\lambda_x$")
    plt.title(f"Rayleigh Wave Halfspace")
    fid = "rayleigh_hs.png"
else:
    plt.ylabel("depth [km]")
    plt.title(f"$\lambda_x$={lambda_x}km; T={T}s")
    fid = f"rayleigh_{int(T)}.png"

plt.tight_layout()
plt.savefig(f"/Users/chow/Documents/academic/teaching/GEOS604_seismo/lecture_screenshots/surface_waves/{fid}")
plt.show()


