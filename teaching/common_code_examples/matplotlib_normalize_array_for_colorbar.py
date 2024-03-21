import numpy as np
import matplotlib.pyplot as plt

# Random sine wave as data
x = np.arange(0, 2 * np.pi)
y = np.sin(x)
z = np.random.rand(len(x))

# This works to plot color
plt.scatter(x, y, c=z, cmap="jet")  
plt.show()

# But what if we're plotting point by point (e.g., looping on a list)
colors = plt.cm.jet(z)

# Or normalize between 0 1
if False:
    norm = plt.Normalize()
    colors = plt.cm.jet(norm(z))

# Plot piece by piece
for x_, y_, c_ in zip(x, y, colors):
    plt.scatter(x_, y_, c=c_)

plt.show()
