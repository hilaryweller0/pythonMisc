#!/usr/bin/python

# Horizontal advection over Schar mountain
from __future__ import division
import pylab as pl
import sys

if len(sys.argv) > 1:
    h_offset = sys.argv[1]
else:
    h_offset = 500

pl.ion()

def h(x):
    x = x + h_offset
    "The mountain height as a function of x"
    h0 = 3e3
    a = 25e3
    lam = 8e3
    return pl.where(abs(x) <= a,
                    h0*pl.cos(pl.pi*x/lam)**2 * pl.cos(pl.pi*x*0.5/a)**2,
                    0)

def sigmaHeights(x,Z,zmax):
    "The height of sigma coordinate Z"
    hx = h(x)
    return hx + Z/zmax*(zmax - hx)

# Parameters of the grid
dx = 1e3
dZ = 500.
dt = 25.
xmin = -150e3
xmax = 150e3
zmin = 0.
zmax = 25e3
tEnd = 10e3
nx = int((xmax - xmin)/dx)
nz = int((zmax - zmin)/dZ)
nt = int(tEnd/dt)

# parameters of the flow
u0 = 10.
z1 = 4e3
z2 = 5e3

# 1d arrays (x and z locations for different variables
Xrho = xmin + (pl.arange(nx)+0.5)*dx
Zrho = zmin + (pl.arange(nz)+0.5)*dZ
xPsi = xmin + pl.arange(nx+1)*dx
ZPsi = zmin + pl.arange(nz+1)*dZ
zPsi = sigmaHeights(pl.reshape(xPsi,[nx+1,1]),pl.reshape(ZPsi,[1,nz+1]),zmax)

# 2d arrays of x and z
xrho = pl.reshape(Xrho,[nx,1])*pl.ones([1,nz])
zrho = sigmaHeights(pl.reshape(Xrho,[nx,1]),pl.reshape(Zrho,[1,nz]),zmax)

# Set up the 2D arrays

# Setup and plot the initial tracer field
rho = pl.zeros([nx, nz])
Ax = 25e3
Az = 3e3
x0 = -50e3
z0 = 9e3
for i in range(nx):
    x = Xrho[i]
    for k in range(nz):
        z = zrho[i,k]
        r = pl.sqrt(((x-x0)/Ax)**2 + ((z-z0)/Az)**2)
        if r <= 1:
            rho[i,k] = pl.cos(pl.pi*r/2)**2

pl.figure(1)
pl.clf()
pl.contour(xrho, zrho, rho, pl.arange(-1.1,1.1,0.1))
pl.colorbar()

# The Jacobian of the coordinate transform for the sigma coordinates
J = zmax/(zmax - pl.reshape(h(Xrho), [nx,1]))*pl.ones([1,nz])

# Stream function at doubly staggered grid points
Psi = pl.zeros([nx+1, nz+1])
for i in range(nx+1):
    x = xPsi[i]
    for k in range(nz+1):
        z = zPsi[i,k]
        if (z <= z1):
            Psi[i,k] = 0
        elif (z <= z2):
            Psi[i,k] = -0.5*u0*(z-z1-(z2-z1)/pl.pi*pl.sin(pl.pi*(z-z1)/(z2-z1)))
        else:
            Psi[i,k] = -0.5*u0*(2*z - z1 - z2)

# Setup leapfrog time-stepping (backwards one for the old old time step)
rhoNew = rho.copy()
rhoOld = rho + 0.5*dt*u0/dx\
         *(pl.roll(rho,-1,axis=0) - pl.roll(rho,1,axis=0))

# Time steps to store for plotting
rhoInit = rho.copy()
rhoMid = rho.copy()

# loop over all times
for it in range(0,nt):
    # loop over space (away from the boundaries, no change on boundaries)
    for i in range(1,nx-1):
        for k in range(1,nz-1):
            rhoNew[i,k] = rhoOld[i,k] - J[i,k]*dt/(dx*dZ)*\
            (\
                - (Psi[i+1,k+1]-Psi[i+1,k])*rho[i+1,k] \
                + (Psi[i,k+1]-Psi[i,k])*rho[i-1,k] \
                + (Psi[i+1,k+1]-Psi[i,k+1])*rho[i,k+1] \
                - (Psi[i+1,k]-Psi[i,k])*rho[i,k-1] \
            )
    # update old arrays
    rhoOld = rho.copy()
    rho = rhoNew.copy()
    
    # store mid value if we have arrived there
    if it*dt == tEnd/2:
        rhoMid = rho.copy()
    
    # plot every few time-steps
    if ((it+1)%10 == 0):
        pl.figure(1)
        pl.clf()
        pl.contour(xrho, zrho, rho, pl.arange(-1.1,1.1,0.1))
        pl.draw()

# Final plot with start, end and mid time-steps
pl.figure(1)
pl.clf()
pl.contour(xrho, zrho, rho,
           pl.concatenate((pl.arange(-2,0,0.1), pl.arange(0.1,2,0.1))),
           colors = 'k')
pl.contour(xrho, zrho, rhoMid,
           pl.concatenate((pl.arange(-2,0,0.1), pl.arange(0.1,2,0.1))),
           colors = 'k', hold='on')
pl.contour(xrho, zrho, rhoInit,
           pl.concatenate((pl.arange(-2,0,0.1), pl.arange(0.1,2,0.1))),
           colors = 'k', hold='on')
pl.draw()
pl.ioff()
pl.show()

