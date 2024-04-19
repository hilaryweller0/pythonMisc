import numpy as np

def advection(phi, fluxFunction, c):
    flux = fluxFunction(phi)
    return -c*(flux - np.roll(flux,1))


def QUICKflux(phi):
    return 0.125*(3*np.roll(phi,-1) + 6*phi - np.roll(phi,1))

def cubicFlux(phi):
    return 1/6*(2*np.roll(phi,-1) + 5*phi - np.roll(phi,1))

def cubicFluxCorr(phi):
    dx = 1/len(phi)
    gradUp = (np.roll(phi,-1) - np.roll(phi,1))/(2*dx)
    gradf = (np.roll(phi,-1) - phi)/dx
    
    print('gradUp =', gradUp)
    print('gradf =', gradf)
    
    return phi + 1/6*dx*(2*gradUp + gradf)

def QUICKfluxCorr(phi):
    dx = 1/len(phi)
    gradUp = (np.roll(phi,-1) - np.roll(phi,1))/(2*dx)
    gradf = (np.roll(phi,-1) - phi)/dx
    
    print('gradUp =', gradUp)
    print('gradf =', gradf)
    
    return phi + 1/4*dx*(gradUp + gradf)

