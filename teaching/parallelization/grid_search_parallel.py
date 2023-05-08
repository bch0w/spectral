import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor

x = np.arange(-1 * np.e, np.e, .001)
y = np.arange(-1 * np.e, np.e, .001)

def gaussian_2d(x_):
    """Calculate 2D Gaussian array for a given X value 
    """
    return np.exp(-1 * x_**2 - y**2)

# Here, the order of returned results matters, so we use
# the map() function
z=[]
with ProcessPoolExecutor() as executor:
    for z_ in executor.map(gaussian_2d, x):
        z.append(z_)

z = np.array(z)
z = z.reshape((len(x), len(y)))
plt.imshow(z)
plt.savefig("2d_gaussian.png")
