import matplotlib.pyplot as plt
import math
import numpy as np


def processFunction(t, tmax=0., sigma=1, C=10., offset=10.):
   o = C * np.exp(-1/2 * (t-tmax)**2 / sigma**2) + offset
   return o

t = np.linspace(0, 10, 50)
f1 = processFunction(t, 2, 0.1)
f2 = processFunction(t, 5, 1)

fig, ax = plt.subplots()
ax.plot(t, f1, c="b", label="f1")
ax.plot(t, f2)
ax.legend()
plt.show()
