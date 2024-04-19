# -*- coding: utf-8 -*-
"""
Mis-matched functions for stability, think FB will have the same stability as Eu (from looking at eqn)
so only Eu is presented.
"""
from __future__ import absolute_import, division
import numpy as np
import matplotlib.pyplot as plt







def eigenSolver(lorenzParams):
    '''
    finds the eigenvales of the linearised lorenz matrix
    '''
    sigma = lorenzParams['sigma']
    rho  = lorenzParams['rho']
    beta = lorenzParams['beta']

    A =np.matrix([[-sigma, sigma,0],[rho, -1,0],[0,0,-beta]])
    e =np.linalg.eigvals(A) 
    return e
    

def numericalAmp(e,dt):
    '''Takes two inputs, e, a set of eigen values and dt a time step.
    It returns the amplification factor corresponding to these eign values.
    '''
    u = dt*e
    #calculate amplification factor
    A_euler = u
    #musRk is the 3 terms that needed to be smaller than mu(denoted u  or A_euler here) for stability
    musRK = 0.5*u**2 +(1/6)*u**3+(1/24)*u**4
    return musRK, A_euler 
    
    

def Main_Plots(lorenzParams,dt_range,scheme):  
    '''
    Going to try and find the stability of the lorenz attractor
    For each scheme there is a different dt_range, execution and plotting protocol
    Poorly organised.
    '''
    B = lorenzParams['beta']
    plt.figure(figsize=(5,4))
    #just input dt as zero as Im not using it in that function for now
    e =eigenSolver(lorenzParams)
    print e
    allmu =[]
    A_E = []
    for dt in dt_range:
        
        if scheme == 'Eu':
            #############Euler#########
            #find max eigenvalue
            #e_abs = [abs(i) for i in e]
            k = max(e)
            plt.plot(dt,k*dt, '.',c ='k')
        
        
        elif scheme =='RK':
            ###########RK4###########
            k = max(e)
            mu=k*dt
            plt.plot(dt,0.5*mu+ (mu**2)/6.+ (mu**3)/24., '.',c='k')
        
   

    
    plt.xlabel('$\Delta t$')
    
    
    
    if scheme == 'Eu':
        plt.ylabel('max(|lambda|) dt') #Euler
        plt.savefig('Eu_FB_stab.pdf')
    
    if scheme == 'RK':
        plt.ylabel('1/2 $ \lambda \Delta$ t +1/6 $ \lambda \Delta$ t$^{2}$ +1/24 $ \lambda \Delta$ t$^{3}$')
        #add horizontal line at y=1 to show limit of numerical stability 
        plt.axhline(y=1, c='b')
        plt.savefig('RK_stab.pdf')
    
    plt.show()

    
def main():

   lorenzParam = {}
   execfile("lorenzParams.py", lorenzParam)
   
   dt_rangeRK = np.linspace(0.0001,0.15,100)#RK4
   dt_rangeEu = np.linspace(0.000001,0.01,80)#Euler
   Main_Plots(lorenzParam, dt_rangeRK, 'RK')
   Main_Plots(lorenzParam, dt_rangeEu, 'Eu')
    
main()
    
    

