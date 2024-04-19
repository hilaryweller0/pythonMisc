import numpy as np
import math
import matplotlib.pyplot as plt

x = np.arange(0,1,0.1)
x = x*2*math.pi
plt.plot(x, np.sin(x))
plt.show()
plt.ion()

plt.clf()
A = np.matrix( [[1,2,3],[11,12,13],[21,22,23]])
C = A.I
B = np.asarray(C)
plt.pcolor(B)
plt.colorbar()

