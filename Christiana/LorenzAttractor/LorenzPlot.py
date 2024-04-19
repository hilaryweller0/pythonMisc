# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 10:53:01 
Contains various plotting functions for solutions and errors of lorenz attractor
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np




def errorPlot(diff, lorenzParam,name,dttype,nt,xLower,xUpper):
    '''
    Plots the error for 1 (constant time-step) case
    '''
    dt = lorenzParam['dt']
  
    duration = dt*nt
    difflog = [np.log(x) for x in diff]
    t = np.linspace(0,duration,len(diff))
    fig, ax1 = plt.subplots()
    ax1.plot(t,diff, 'k')
    ax1.set_xlabel('duration, t')
    ax1.set_ylabel('RMS error in position:'+str(name)+' and RK4 High Res')
    
    ax2 = ax1.twinx()
    ax2.plot(t,difflog, c='b')
    ax2.set_ylabel('log of RMS error in position:'+str(name)+' and RK4 High Res', color='k')
    for tl in ax2.get_yticklabels():
        tl.set_color('b')
        
    plt.axvline(x = xLower, ls = '--',color = 'lightslategrey')
    plt.axvline(x = xUpper, ls = '--', color = 'lightslategrey')
    plt.savefig('LongRunPlot_' + str(name) + '.pdf')
    plt.show()
    




def errorCompar(d_const,d_var,dt,nt,name,dt_var,var_type):
    '''
    plots two errors one constant one variable
    d_const is the errors of the constant case, d_var - erors of variable case
    dt_var is the variable case time steps
    '''
    
    duration = dt*nt
    #const case time axis
    t = np.linspace(0,duration,len(d_const))
    t_var = np.cumsum(dt_var)

    plt.plot(t,d_const, 'k')
    plt.plot(t_var,d_var,'b')
    plt.xlabel('duration, t')
    plt.ylabel('RMS error in position wtr RK: '+str(name)+' and '+str(name)+ 'variable')

    
    plt.savefig('ComparisionPlot_'+str(name)+'variable_'+str(var_type)+'.pdf')
    plt.show()
    
    
    
    
    


def lorentzPlotting(x,y,z,lorenzParam,nt):
    '''
    Plots the solution of the lorenz equations
    '''

    sigma = lorenzParam['sigma']
    rho = lorenzParam['rho']
    beta = lorenzParam['beta']
    dt = lorenzParam['dt']
    nt = len(x)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x,y,z)
    plt.gca().patch.set_facecolor('white')
    
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    fig.savefig('Attractor_' + str(nt) + 'timesteps.pdf')
    fig.show()
    
    
def ICplot(x,y,z,lorenzParam,xs,ys,zs):
    '''
    plots the chosen 100 initial contions x,y,z out of the larger section xs,ys,zs
    xyz are points on a xs ys zs solution
    '''
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xs,ys,zs,c = 'royalblue',alpha = 0.6)
    ax.scatter(x,y,z,'o', alpha =1.0)
    
    plt.gca().patch.set_facecolor('white')
    
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')  
    fig.savefig('IC_distribution_plot.pdf')
    fig.show()

