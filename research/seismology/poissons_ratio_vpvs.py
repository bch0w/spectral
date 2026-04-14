
"""
Plot Vp/Vs ratio against Poisson's Ratio for reference
"""
import numpy as np
import matplotlib.pyplot as plt


vpvs = np.linspace(1.3, 2.3, 100)
poissons = 0.5 * (vpvs ** 2 - 2) / (vpvs ** 2 - 1)

plt.plot(vpvs, poissons, "k-", linewidth=2)
plt.axvline(1.73, color="r")
plt.xlabel("Vp/Vs")
plt.ylabel("Poisson's Ratio")
plt.title("Vp/Vs Ratio against Poisson's Ratio\n"
          "(red line = Poisson's solid)")
plt.grid()
plt.show()
