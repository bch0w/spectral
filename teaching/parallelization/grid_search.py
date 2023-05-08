import numpy as np
import matplotlib.pyplot as plt
z = []
x = np.arange(-1 * np.e, np.e, .001)
y = np.arange(-1 * np.e, np.e, .001)
for x_ in x:
    for y_ in y:
        z_ = np.exp(-1 * x_**2 - y_**2)  # 2D Gaussian
        z.append(z_)

z = np.array(z)
z = z.reshape((len(x), len(y)))
plt.imshow(z)
plt.savefig("test.png")
