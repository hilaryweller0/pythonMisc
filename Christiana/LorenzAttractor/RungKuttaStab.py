# -*- coding: utf-8 -*-
"""
Created on Sun May  8 11:48:31 2016

@author: christiana
"""
"""
Going to try and plot the stability region of the 4th order Runge Kutta Scheme
"""

import numpy as np
import matplotlib.pyplot as plt

gamma_dt = np.linspace(-3,3,40)
omega_dt = np.linspace(-3,3,40)
#create the coord space
[gam,omeg] = np.meshgrid(gamma_dt,omega_dt)
#create complex value at each point
z = gam +1j*omeg
# find the amplification factor at each point
A = 1 + z +0.5*z**2 + (z**3)/6. +(z**4)/24.
mag_A = abs(A)
plt.contour(gamma_dt,omega_dt,mag_A,[1])

plt.axvline(x=0, ymin=0, ymax = 3,color='k')
plt.xlabel('gamma dt')
plt.ylabel('omega dt')
