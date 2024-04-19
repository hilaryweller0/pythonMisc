import numpy as np
import matplotlib.pyplot as plt
from advectionFluxes import *
from initialConditions import *
from RKschemes import *

# Problem parameters
nx = 12

# Derived parameters
dx = 1/nx

x = np.arange(0.5*dx,1+0.5*dx,dx)

# Initial conditions
T = x**3
divTExact = 3*x**2

print('Cubic flux = ', cubicFluxCorr(T))

print('QUICK flux = ', QUICKfluxCorr(T))

divT = -advection(T, cubicFluxCorr, 1/dx)
print('Div error ', divT - divTExact)

divT = -advection(T, QUICKfluxCorr, 1/dx)
print('Div error ', divT - divTExact)

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

