# Various useful diagnostics for solving the Diffusion equation
from __future__ import division
import pylab as pl
from scipy import special

def analyticErf(x, Kt, alpha, beta):
    "The analytic solution of the 1d diffusion equation with diffuions"
    "coeffienct K at time t assuming top-hat initial conditions which are"
    "one between alpha and beta and zero elsewhere"
    
    phi = 0.5 * special.erf((x-alpha)/pl.sqrt(4*Kt))  \
        - 0.5 * special.erf((x-beta )/pl.sqrt(4*Kt))
    return phi

def fourierPower(signal):
    "Calculates the Fourier power spectrum of the given array"
    
    A = pl.fft(signal)
    n = (len(signal)+1)//2
    power = pl.zeros(n)
    for i in range(n):
        power[i] = 4*(A[i+1].real**2 + A[i+1].imag**2)
    return power

def l2ErrorNorm(phi, phiExact):
    "Calculates the l2 error norm (RMS error) of phi in comparison to"
    "phiExact"
    
    #remove one of the end points
    phi = phi[0:-1]
    phiExact = phiExact[0:-1]
    
    # calculate the error and the error norms
    phiError = phi - phiExact
    l2 = pl.sqrt(sum(phiError**2)/sum(phiExact**2))

    return l2

