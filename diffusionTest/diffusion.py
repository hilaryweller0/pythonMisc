#!/usr/bin/python

# Outer code for setting up the diffusion problem and calling the
# function to diffusion.
from __future__ import division
import pylab as pl

# all the diffusion schemes and initial conditions
execfile("diffusionSchemes.py")
execfile("initialConditions.py")
execfile("diagnostics.py")

def diffusion(xmin = 0, xmax = 1, nx = 40, nt = 40, K = 1e-3, dt = 0.1):
    "Diffuse between x = xmin and x = xmax split over 40 spatial steps"
    "with diffusion coeffient K, time step dt for nt time steps"

    # Derived parameters
    dx = (xmax - xmin)/nx
    Kdtbydx2 = K*dt/dx**2   # Non-dimensional diffusion coeffient
    print "non-dimensional diffusion coefficient = ", Kdtbydx2
    print "dx = ", dx, " dt = ", dt, " nt = ", nt
    print "end time = ", nt*dt

    # spatial points for plotting and for defining initial conditions
    x = pl.linspace(xmin, xmax, nx+1)

    # initial conditions
    phiOld = topHat(x, 0.4, 0.6)
    # Analytic solution (of top-hat profile)
    phiAnalytic = analyticErf(x, K*dt*nt, 0.4, 0.6)

    # Diffusion using FTCS and BTCS
    phiFTCS = FTCS(phiOld.copy(), Kdtbydx2, nt)
    phiBTCS = BTCS(phiOld.copy(), Kdtbydx2, nt)

    # plot options
    font = {'size'   : 20}
    pl.rc('font', **font)
    #pl.ion()

    # plot the results in comparison to the initial conditions and analytic
    pl.figure(1)
    pl.clf()
    #pl.plot(x, phiOld,     'k--', label='Initial conditions')
    pl.plot(x, phiAnalytic,'k',label='Exact solution')
    pl.plot(x, phiFTCS,     'r', label='FTCS')
    pl.plot(x, phiBTCS,     'b', label='BTCS')
    pl.legend(loc=4)
    pl.xlabel('x')
    #pl.ylabel('$\phi$')

    pl.figure(2)
    pl.clf()
    pl.plot(x, phiFTCS-phiAnalytic, 'r', label='FTCS')
    pl.plot(x, phiBTCS-phiAnalytic,     'b', label='BTCS')
    pl.legend(loc='best')
    pl.xlabel('x')
    #pl.ylabel('$\phi$ error')

    print "FTCS l2 error norm = ", l2ErrorNorm(phiFTCS, phiAnalytic)
    print "BTCS l2 error norm = ", l2ErrorNorm(phiBTCS, phiAnalytic)

    # plot options
    font = {'size'   : 16}
    pl.rc('font', **font)

    # power spectra
    powerInit = fourierPower(phiOld[0:-1])
    powerFinal = fourierPower(phiAnalytic[0:-1])
    powerFTCS = fourierPower(phiFTCS[0:-1])
    powerBTCS = fourierPower(phiBTCS[0:-1])
    kdx = pl.arange(dx, dx*len(powerInit), dx)
    pl.figure(3)
    pl.clf()
    pl.semilogy(kdx,powerInit[1:], 'k--', label='Initial conditions')
    pl.semilogy(kdx,powerFinal[1:],'k',   label='Exact solution')
    pl.semilogy(kdx,powerFTCS[1:], 'r',   label='FTCS')
    pl.semilogy(kdx,powerBTCS[1:], 'b',   label='BTCS')
    pl.legend(loc='best')
    pl.xlabel('$k\Delta x$')
    pl.ylabel('power')

diffusion(nx = 40, nt = 40, dt = 0.1)
pl.figure(1)
pl.savefig('phi40.pdf')
pl.figure(2)
pl.savefig('phiError40.pdf')
pl.figure(3)
pl.savefig('phiPower40.pdf')
