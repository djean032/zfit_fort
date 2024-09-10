import numpy as np
import matplotlib.pyplot as plt

x, y = np.loadtxt('intensities.dat', usecols=(0, 1), unpack=True)
plt.plot(x, y)
plt.show()
