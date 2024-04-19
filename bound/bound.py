import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-1.,1.,101)
phi = x
eps = 0.1

bound1 = np.where(phi < eps, eps, phi)
phiBar = np.mean(bound1)
pos0 = np.where(-phi < 0, 0, -phi)
bound2 = np.where(phi < phiBar*pos0,  phiBar*pos0, phi)
bound = np.where(bound2 < eps, eps, bound2)

plt.plot(x, phi, label='phi')
plt.plot(x, bound, label='bounded')
plt.plot(x, eps*np.ones_like(x), label='bound')
plt.axvline(eps)
plt.axvline(-eps)
plt.legend()
plt.show()

