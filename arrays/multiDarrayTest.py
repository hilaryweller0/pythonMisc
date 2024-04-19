import numpy as np
import matplotlib.pyplot as plt

nx = 41
x = np.linspace(0,2*np.pi, nx)

sinx = np.zeros([2,nx,3])
sinx[0,:,0] = np.sin(x)
sinx[0,:,1] = np.sin(2*x)
sinx[0,:,2] = np.sin(3*x)
sinx[1,:,0] = np.sin(x/2)
sinx[1,:,1] = np.sin(x/3)
sinx[1,:,2] = np.sin(x/4)

plt.plot(x, sinx[0,:,0])
plt.plot(x, sinx[0,:,1])
plt.plot(x, sinx[0,:,2])
plt.plot(x, sinx[1,:,0])
plt.plot(x, sinx[1,:,1])
plt.plot(x, sinx[1,:,2])
plt.show()

