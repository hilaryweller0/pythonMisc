import numpy as np

def RK3(y, dydt, flux, dt):

    F = dydt(y, flux, dt)
    yNew = y + dt*F
    
    F += dydt(yNew, flux, dt)
    yNew = y + 0.25*dt*F
    
    F += 4*dydt(yNew, flux, dt)
    yNew = y + dt/6*F
    
    return yNew


def RK_compact(y, dydt, flux, dt, A):
    F = A[0][0]*dydt(y, flux, dt)
    yNew = y + A[0][1]*F
    for it in range(1, len(A)):
        F += A[it][0]*dydt(yNew, flux, dt)
        yNew = y + A[it][1]*F
    return yNew
