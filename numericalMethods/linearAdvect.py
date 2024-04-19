import numpy as np
import matplotlib.pyplot as plt
from advectionFluxes import *
from initialConditions import *
from RKschemes import *

# Problem parameters
nx = 40
nt = 10
c = 0.01
u = 1
RKscheme = [[1,1], [1,.25], [4,1/6]]

# Derived parameters
dx = 1/nx
dt = c*dx/u
Tmax = nt*dt

print('c =', c, 'dx =', dx, 'dt =', dt, 'Tmax =', Tmax)

x = np.arange(0,1,dx)

# Initial conditions
initial = cosWave
phi = initial(x)
# Analytic solution
phiExact = initial(x - u*Tmax)

# Advection over time
for it in range(nt):
    #FE phi = phi + advection(phi, QUICKflux, c)
    phi = RK_compact(phi, advection, cubicFlux, c, RKscheme)

plt.plot(x, phi, label='cubic, RK3')
plt.plot(x, phiExact, label='Analytic')
plt.legend()
plt.show()

# Convergence with resolution
nxs = [20, 40, 80]
nts = [5, 10, 20]
dxs = np.zeros(len(nxs))
errors = np.zeros(len(nxs))
for i in range(len(nxs)):
    nx = nxs[i]
    dx = 1/nx
    dxs[i] = dx
    dt = c*dx/u
    nt = nts[i]
    Tmax = nt*dt
    x = np.arange(0, 1, dx)
    phi = initial(x)
    phiExact = initial(x - u*Tmax)
    
    # Advection over time
    for it in range(nt):
        phi = RK_compact(phi, advection, cubicFlux, c, RKscheme)
    errors[i] = l2error(phi, phiExact)

print(errors)

