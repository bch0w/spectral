"""
Plot tangent without including the drops along the asymptotes
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, tan


x = np.linspace(- pi / 2, 4 * pi, 1000)

for i, n in enumerate([1, 0.2, 0.1]):
    y = tan(n * x)
    y[:-1][np.diff(y) < 0] = np.nan   # remove negative gradients
    plt.plot(x, y, c=f"C{i}", label=f"n={n}")
plt.ylim([-6, 6])
plt.xlim([x.min(), x.max()])
plt.xlabel("x")
plt.ylabel("tan(nx)")
plt.legend()
plt.grid()
plt.show()

