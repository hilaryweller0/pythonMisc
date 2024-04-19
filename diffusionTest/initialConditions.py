# Initial conditions for the diffusion equation
from __future__ import division
import pylab as pl
import sys

def topHat(x,alpha,beta):
    "Function defining a top hat as a function of position, x"
    "chooses 1 where condition is true, else chooses zeros"
    "very close to the limits it choses 0.5"
    
    # smallest float
    eps = sys.float_info.epsilon
    
    phi = pl.where(
            (x > alpha+eps) & (x < beta-eps), 1., 
          pl.where(
            pl.absolute(x-beta)<eps, 0.5, 
          pl.where(
            pl.absolute(x-alpha)<eps, 0.5,
            0)))
    return phi

