# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 11:14:18 2016

@author: christiana
"""
import numpy as np
import matplotlib.pyplot as plt

#Lorenz Parameters for a 'chaotic solution'
beta  = 8./3.
rho   = 28.
sigma = 10.
#constant time step 
dt=0.002

#the number of timesteps if this length required to access the region of exponential error growth
#is different for Euler and RK4 and BF so the require different nts
nt_RK = 6500#5000 #7000 timestep <0.07 for dual attractor, <0.2 for solution 
nt_FB =1400#1500 #3500 timestep < 0.07, works for 0.06
nt_Eu = 500#350 #2500 350 timestep = <0.03, works for 0.02

#variable time step making sure duration and number of timesteps is the same
#as the constant timestep case

############################## RK 4 #######################################


#variables for the *a^n case
#Number from solving via wolfram, added 5 on end to make both durations match
variable_dtRK = [(0.00002)*1.000996679189304**i for i in range(nt_RK)]
variable_timeRK = np.cumsum(variable_dtRK)
const_timeRK = np.cumsum([dt]*nt_RK)
var_typeRK = 'multiplication'



'''
#variables for the _n*a case
variable_dtRK = [(0.00008)+i*(5.90860132328e-7)for i in range(nt_RK)]
variable_timeRK = np.cumsum(variable_dtRK)
const_timeRK = np.cumsum([dt]*nt_RK)
var_typeRK = 'addition'
'''

print variable_timeRK[-1], const_timeRK[-1]
########################## Euler Forward ###################################

#variables for the *a^n case
variable_dtEu = [(0.00002)*1.01304892418836**i for i in range(nt_Eu)]
variable_timeEu = np.cumsum(variable_dtEu)
const_timeEu = np.cumsum([dt]*nt_Eu)
var_typeEu = 'multiplication'


'''
#variables for the +n*a case
variable_dtEu = [(0.00008)+i*(7.69539078156e-6) for i in range(nt_Eu)]
variable_timeEu = np.cumsum(variable_dtEu)
const_timeEu = np.cumsum([dt]*nt_Eu)
var_typeEu = 'addition'

'''


print variable_timeEu[-1], const_timeEu[-1]

######################### Backwards forwards ###############################

#variables for the *a^n case
variable_dtFB = [(0.00002)*1.004637385427204**i for i in range(nt_FB)]
variable_timeFB = np.cumsum(variable_dtFB)
const_timeFB = np.cumsum([dt]*nt_FB)
var_typeFB = 'multiplication'



'''

#variables for the +n*a case
variable_dtFB = [0.00008+i*(2.74481772695e-6) for i in range(nt_FB)]
variable_timeFB = np.cumsum(variable_dtFB)
const_timeFB = np.cumsum([dt]*nt_FB)
var_typeFB = 'addition'

'''
print variable_timeFB[-1], const_timeFB[-1]
#variable_dt = [(4*dt/10)+n*2*dt/100000 for n in range(0,nt)], need nt=40001
'''
variable_dt1 = [0.0002]*(nt/2) 
variable_dt2 = [0.002]*(nt/2)
variable_dt = variable_dt1+variable_dt2



t = np.cumsum(variable_dt)
th = [dt/10]*(nt*10)
'''