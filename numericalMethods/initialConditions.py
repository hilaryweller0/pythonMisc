import numpy as np

def cosBell(x, a = 0., b = 0.5):
    return np.where((x>a) & (x<b), 1-np.cos(2*np.pi*(x-a)/(b-a)), 0.)

def cosWave(x, k=1):
    return np.cos(2*np.pi*x*k)

def l2error(phi, phiExact):
    return np.sqrt(sum((phi - phiExact)**2)/sum(phi**2))

