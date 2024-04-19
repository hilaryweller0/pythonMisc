from __future__ import division
import pylab as pl
import scipy.linalg

def FTCS(phiOld, k, nt):
    "Diffusion of profile in phiOld using FTCS using non-dimensional"
    "diffusion coefficient, k"
    
    # new time-step array for phi
    phi = pl.zeros_like(phiOld)
    

    return phi


def BTCS(phi, k, nt):
    "Diffusion of profile in phi using BTCS using non-dimensional"
    "diffusion coefficient, k"
    
    N = len(phi)

    return phi


